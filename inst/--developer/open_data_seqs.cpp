#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <algorithm>

using dbl_vec = std::vector<double>;
using str_vec = std::vector<std::string>;

class CSVReader
{
public:
    std::vector<dbl_vec> data_seqs;
    std::unordered_map<std::string, size_t> idToIndex;
    str_vec orderedIds;

    // Read CSV file and allocate memory for rows
    bool readCSV(const std::string &filename)
    {
        data_seqs.clear();
        idToIndex.clear();
        orderedIds.clear();

        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        std::string line;
        size_t rowIndex = 0;

        while (std::getline(file, line))
        {
            // Skip empty lines
            if (line.empty())
                continue;

            // Parse the line
            std::stringstream ss(line);
            std::string token;
            bool firstToken = true;
            std::string id;
            std::vector<double> values;

            while (std::getline(ss, token, ','))
            {
                // Remove leading/trailing whitespace
                token.erase(0, token.find_first_not_of(" \t\r\n"));
                token.erase(token.find_last_not_of(" \t\r\n") + 1);

                if (firstToken)
                {
                    // First token is the ID
                    id = token;
                    firstToken = false;
                }
                else
                {
                    // Subsequent tokens are numerical values
                    try
                    {
                        double value = std::stod(token);
                        values.push_back(value);
                    }
                    catch (const std::exception &e)
                    {
                        std::cerr << "Warning: Could not parse value '" << token
                                  << "' in row " << rowIndex + 1 << std::endl;
                    }
                }
            }

            // Store the data_seqs
            if (!id.empty())
            {
                // Check for duplicate IDs
                if (idToIndex.find(id) != idToIndex.end())
                {
                    std::cerr << "Warning: Duplicate ID '" << id
                              << "' found. Overwriting previous entry." << std::endl;
                    size_t existingIndex = idToIndex[id];
                    data_seqs[existingIndex] = std::move(values);
                }
                else
                {
                    // Add new entry
                    idToIndex[id] = rowIndex;
                    orderedIds.push_back(id);
                    data_seqs.push_back(std::move(values));
                    rowIndex++;
                }
            }
        }

        file.close();

        return true;
    }

    // Get values by row index (returns nullptr if not found)
    const dbl_vec *getRowByIndex(size_t index) const
    {
        if (index < data_seqs.size())
            return &data_seqs[index];
        return nullptr;
    }

    // Get values by ID (returns nullptr if not found)
    const dbl_vec *getRowById(const std::string &id) const
    {
        auto it = idToIndex.find(id);
        if (it != idToIndex.end())
            return &data_seqs[it->second];
        return nullptr;
    }

    // Get index for a given ID (returns -1 if not found)
    size_t getIndexById(const std::string &id) const
    {
        auto it = idToIndex.find(id);
        if (it != idToIndex.end())
            return static_cast<size_t>(it->second);
        return -1;
    }

    // Get ID for a given index (returns empty string if not found)
    std::string getIdByIndex(size_t index) const
    {
        if (index < orderedIds.size())
            return orderedIds[index];
        return "";
    }

    // Print all data_seqs (for debugging)
    void printData() const
    {
        for (size_t i = 0; i < orderedIds.size(); ++i)
        {
            std::cout << "ID: " << orderedIds[i] << ", Values: ";
            for (const auto &val : data_seqs[i])
            {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }            

    // Iterator support for range-based loops
    class Iterator
    {
    private:
        const CSVReader *reader;
        size_t index;

    public:
        Iterator(const CSVReader *r, size_t i) : reader(r), index(i) {}

        bool operator!=(const Iterator &other) const
        {
            return index != other.index;
        }

        Iterator &operator++()
        {
            ++index;
            return *this;
        }

        std::pair<std::string, const std::vector<double> *> operator*() const
        {
            return std::make_pair(reader->orderedIds[index], &reader->data_seqs[index]);
        }
    };

    Iterator begin() const { return Iterator(this, 0); }
    Iterator end() const { return Iterator(this, data_seqs.size()); }
};

// Example usage
int main()
{
    CSVReader reader;

    // Create sample CSV file for testing
    std::ofstream testFile("test_data.csv");
    testFile << "id_1,3.1,7.2,3.2\n";
    testFile << "id_2,1.2,1.1,5.0,4.1,1.2\n";
    testFile << "alpha,2.5,3.7\n";
    testFile << "42,9.9,8.8,7.7,6.6\n";
    testFile.close();

    if (reader.readCSV("test_data.csv"))
    {
        std::cout << "\n=== Full Data ===" << std::endl;
        reader.printData();

        // Access by ID
        std::cout << "\n=== Access by ID ===" << std::endl;
        const auto *row = reader.getRowById("id_2");
        if (row)
        {
            std::cout << "Values for 'id_2': ";
            for (double val : *row)
            {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }

        // Access by index
        std::cout << "\n=== Access by Index ===" << std::endl;
        std::string id = reader.getIdByIndex(2);
        if (!id.empty())
        {
            std::cout << "ID at index 2: " << id << std::endl;
            const auto *values = reader.getRowByIndex(2);
            if (values)
            {
                std::cout << "Values: ";
                for (double val : *values)
                {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
        }

        // Range-based loop (C++14 compatible)
        std::cout << "\n=== Range-based Loop ===" << std::endl;
        for (const auto &pair : reader)
        {
            std::cout << "ID: " << pair.first << ", first value: " << pair.second->at(0) << std::endl;
        }
    }

    return 0;
}