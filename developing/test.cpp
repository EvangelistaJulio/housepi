#include <iostream>
#include <vector>
#include <unordered_map>

std::vector<double> globalVector = {1.1, 2.2, 3.3, 4.4}; // Global vector of double

int main() {

    // Define an unordered map with keys as strings and values as pointers to vector entries
    std::unordered_map<std::string, double*> myMap;

    // Populate the map with pointers to vector entries
    myMap["key1"] = &globalVector[0];
    myMap["key2"] = &globalVector[1];
    myMap["key3"] = &globalVector[2];
    myMap["key4"] = &globalVector[3];

    // Print the initial value associated with key 'key2'
    std::cout << "Initial value associated with key 'key2': " << *(myMap["key2"]) << std::endl;

    // Modify the global vector
    globalVector[1] = 5.5;

    // Print the updated value associated with key 'key2'
    std::cout << "Updated value associated with key 'key2': " << *(myMap["key2"]) << std::endl;

    return 0;
}
