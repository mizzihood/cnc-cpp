#include <filesystem>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <ctime>

#include "line.h"

namespace Cnc
{

double parseCoordinateValue(const std::string& token)
{
    return std::stod(token.substr(1));
}

void splitLineComment(std::string& line, std::string& comment)
{
    std::size_t semicolonPosition = line.find(';');
    if (semicolonPosition != std::string::npos)
    {
        // together with the ;
        comment = line.substr(semicolonPosition);
        line.resize(semicolonPosition);
    }
}

void printPosition(std::ostream& os, TVector position)
{
    os << std::setprecision(3) << std::fixed  
       << "X" << position[0] << " " 
       << "Y" << position[1] << " "
       << "Z" << position[2] << " ";
}

int isCoordinate(const std::string& token)
{
    switch (token[0])
    {
    case 'X' : return 0;
    case 'Y' : return 1;
    case 'Z' : return 2;
    default  : return -1;
    }
}

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
    {
        return ""; // no content
    }
    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

void CncCorrection::calculateMatrix(double angleAB, double angleAC, double angleBC)
{
    double x3 = std::cos(angleAC);
    double y3 = (std::cos(angleBC) - x3 * std::cos(angleAB)) / std::sin(angleAB);
    
    TMatrix A;

    A[0][0] = 1; 
    A[0][1] = std::cos(angleAB);
    A[0][2] = std::cos(angleAC);
    
    A[1][0] = 0; 
    A[1][1] = std::sin(angleAB);
    A[1][2] = y3;

    A[2][0] = 0; 
    A[2][1] = 0;
    A[2][2] = std::sqrt(1 - x3 * x3  - y3 * y3);

    double determinant =    +A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
                            -A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                            +A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    double invdet = 1 / determinant;
    m_matrix[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) * invdet;
    m_matrix[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * invdet;
    m_matrix[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * invdet;
    m_matrix[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * invdet;
    m_matrix[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * invdet;
    m_matrix[1][2] = (A[1][0] * A[0][2] - A[0][0] * A[1][2]) * invdet;
    m_matrix[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) * invdet;
    m_matrix[2][1] = (A[2][0] * A[0][1] - A[0][0] * A[2][1]) * invdet;
    m_matrix[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) * invdet;  

    //T = np.matrix([
    //    [1, std::cos(angleAB), std::cos(angleAC)],
    //    [0, std::sin(angleAB), y3],
    //    [0, 0, std::sqrt(1 - x3 * x3  - y3 * y3)]
    //    ])

}

TVector CncCorrection::transform(const TVector& currentPosition)
{
    TVector newPosition;

    for (int i = 0; i < 3; ++i)
    {
        newPosition[i] = 0;
        for (int j = 0; j < 3; ++j)
        {
            newPosition[i] += m_matrix[i][j] * currentPosition[j];
        }
    }

    return newPosition;
}
CncCorrection::CncCorrection()
{
    m_settings["side_x"] = 10;
    m_settings["side_y"] = 10;
    m_settings["diagonal_from_x0y0"] = 14.142135623730951;
    m_settings["height"] = 100;
    m_settings["z_to_x"] = 10;
    m_settings["z_to_y"] = 20;
    m_currentPosition[0] = 0.;
    m_currentPosition[1] = 0.;
    m_currentPosition[2] = 0.;
    m_currentPositionOut[0] = 0.;
    m_currentPositionOut[1] = 0.;
    m_currentPositionOut[2] = 0.;
}

void CncCorrection::loadSettings(std::istream& is)
{
    std::string line;
    while (std::getline(is, line))
    {
        if (line.find('#') != std::string::npos)
        {
            continue;
        }
        
        std::stringstream sLine(line);
        std::string token; 
        // extract the key 
        std::getline(sLine, token, ':');
        std::string key = trim(token);
        // extract the data 
        std::getline(sLine, token);
        std::string data = trim(token);
        auto settingIterator = m_settings.find(key);
        if (settingIterator != m_settings.end())
        {
            settingIterator->second = std::stod(data);
        }
    }
    // solve the parallelogram
    double a = m_settings["side_x"];
    double b = m_settings["side_y"];
    double q = m_settings["diagonal_from_x0y0"];
    // w is addition to a so the b,a+w|q form a right triangle
    double w = (q * q - (a * a + b * b)) / (2 * a);
    double angle_ab = std::acos(w / b);
    
    // calculate angle_ac and angle_bc
    double h = m_settings["height"];
    double x = m_settings["z_to_x"];
    double angle_ac = std::atan2(h, x);
    // angle_bc
    double y = m_settings["z_to_y"];
    double angle_bc = atan2(h, y);
    
    // calculate matrix
    calculateMatrix(angle_ab, angle_ac, angle_bc);
}

void CncCorrection::printSettings(std::ostream& os)
{
    for (const auto& element : m_settings)
    {
        os << element.first << ": " << element.second << std::endl;
    }
}


void CncCorrection::parseLine(std::string& line, std::ostream& os)
{
    std::string comment;
    splitLineComment(line, comment);

    std::stringstream sLine(line);

    std::string token;
    std::vector<std::string> tokens;
    int firstCoordinateTokenIndex = -1;
    int tokenIndex = 0;
    
    while (sLine >> token)
    {
        int coordinateIndex = isCoordinate(token);
        if ( coordinateIndex >= 0)
        {
            m_currentPosition[coordinateIndex] = parseCoordinateValue(token);
            if (firstCoordinateTokenIndex < 0)
            {
                firstCoordinateTokenIndex = tokenIndex;
            }
        }
        else
        {
            tokens.push_back(token);
        }
        ++tokenIndex;
    }
    
    if (comment.size() > 0)
    {
        tokens.push_back(comment);
    }

    TVector positionOut = transform(m_currentPosition);

    tokenIndex = 0;
    for (auto& token : tokens)
    {
        os << token << ' ';
        ++tokenIndex;
        if (tokenIndex == firstCoordinateTokenIndex)
        {
            printPosition(os, positionOut);
        }
    }
    os << std::endl;
}

void CncCorrection::parseFile(std::istream& is, std::ostream& os)
{
    std::string line;

    while (std::getline(is, line))
    {
        parseLine(line, os);
    }
}
} // namespace 

//    os << "";Izvorna datoteka   : " + file_in_path  + "\n")
//    os << "";Ustvarjena datoteka: " + file_out_path + "\n")
//    os << "";Datum obdelave     : " + dt_string     + "\n")


int main(int argc, char* argv[])
{
    std::ifstream fSettings(argv[1], std::ios_base::in);
    Cnc::CncCorrection cnc;
    cnc.loadSettings(fSettings);
    cnc.printSettings(std::cout);
    
    // std::string fileName("..\\mcode\\6_obrez_flat25.nc");
    for (int i = 2; i < argc; ++i)
    {
        const std::filesystem::path fileNameIn = argv[i];
        std::filesystem::path fileNameOut = fileNameIn.parent_path();
        fileNameOut /= "output";
        std::filesystem::create_directories(fileNameOut);
        fileNameOut /= fileNameIn.filename();
        std::cout << "fileNameIn: " << fileNameIn << std::endl 
                  << "fileNameOut: " << fileNameOut << std::endl
                  << "---------------------------------------" << std::endl; 

        std::ifstream mCodeIn(fileNameIn, std::ios_base::in);
        if (!mCodeIn.is_open())
        {
            throw std::invalid_argument("Cannot open file in " + fileNameIn.u8string());
        }
        std::ofstream mCodeOut(fileNameOut, std::ios_base::out); 
        if (!mCodeOut.is_open())
        {
            throw std::invalid_argument("Cannot open file out " + fileNameOut.u8string());
        }

        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        mCodeOut << ";Izvorna datoteka   : " << fileNameIn.u8string() << "\n" 
                 << ";Ustvarjena datoteka: " << fileNameOut.u8string() << "\n" 
                 << ";Datum obdelave     : " << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << "\n";
        cnc.parseFile(mCodeIn, mCodeOut);
    }
    return 0;
}

