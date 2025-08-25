#include "interface.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

/**
 * @brief Reads a CSV file and returns its data as column vectors of doubles.
 * * This function handles CSV files where each value is a double. It reads the
 * entire file into memory, transposes it, and returns it as a vector of columns.
 * It assumes the CSV has a consistent number of columns in each row.
 * * @param filename The path to the CSV file.
 * @return std::vector<std::vector<double>> A vector where each inner vector represents a column of data.
 * @throws std::runtime_error if the file cannot be opened or is empty.
 */
std::vector<std::vector<double>> read_csv_columns(const std::string& filename) {
    // --- 1. Read data row by row ---
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<std::vector<double>> data_rows;
    std::string line;

    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\n\v\f\r") == std::string::npos) {
            continue;
        }

        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            try {
                // Convert string cell to double
                row.push_back(std::stod(cell));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Warning: Could not convert string to double: " << cell << ". Skipping." << std::endl;
                // Optionally, you could push a NaN or a default value
                // row.push_back(std::numeric_limits<double>::quiet_NaN());
            }
        }
        data_rows.push_back(row);
    }
    file.close();

    if (data_rows.empty()) {
        throw std::runtime_error("CSV file is empty or contains no valid data.");
    }

    // --- 2. Transpose from rows to columns ---
    // The number of columns is the size of the first row.
    // We assume all rows have the same number of columns.
    size_t num_cols = data_rows[0].size();
    size_t num_rows = data_rows.size();

    std::vector<std::vector<double>> data_cols(num_cols);

    for (size_t col_idx = 0; col_idx < num_cols; ++col_idx) {
        // Pre-allocate memory for efficiency
        data_cols[col_idx].reserve(num_rows);
        for (size_t row_idx = 0; row_idx < num_rows; ++row_idx) {
            // Basic safety check for ragged CSVs
            if (col_idx < data_rows[row_idx].size()) {
                data_cols[col_idx].push_back(data_rows[row_idx][col_idx]);
            } else {
                // Handle rows with fewer columns than the first one
                // For example, push a default value like 0 or NaN
                data_cols[col_idx].push_back(0.0); 
            }
        }
    }

    return data_cols;
}

int main() {
    auto datacols = read_csv_columns("eos_validated.csv");
    std::vector<double> p = datacols[1];
    std::vector<double> e = datacols[0];

    int n = p.size();
    if (n != e.size()) {
        std::cerr << "Error: The number of pressure and energy density values must match." << std::endl;
        return 1;
    }

    double out[2];
    double eps_c = 100.0;
    double eps_start = 50.0;
    double eps_end = 1200.0;
    int n_stars = 1000;
    double epsc_sol[n_stars];
    double M_sol[n_stars];
    double R_sol[n_stars];

    for (int i = 0; i < 100; i++) {
        qlimr_getMR(p.data(), e.data(), n, eps_c, out);
        // qlimr_getMRdiagram(p.data(), e.data(), n, eps_start, eps_end, n_stars, epsc_sol, M_sol, R_sol);
    }

    return 0;
}