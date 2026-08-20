#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <stdexcept>

struct argument {
    std::string name, short_flag, long_flag, help, defualt_value, data_type;
    bool required;
    std::vector<std::string> values;

    argument(std::string&& name,
             std::string&& short_flag,
             std::string&& long_flag,
             std::string&& help,
             const bool& required,
             std::string&& defualt_value,
             std::string&& data_type);
};

class ArgumentParser {
public:
    static int PRINT_HELP;
    static int ERROR_PARSER;
    static int ERROR_PARSER_INPUTS;
    static int ERROR_PARSER_REQUIRED;
    static int ERROR_PARSER_UNKNOWN;
    static int ERROR_PARSER_INVARG;
    static int ERROR_PARSER_OUTRANGE;
    static int ERROR_PARSER_BOOL;
    static int ERROR_PARSER_CHAR;

    std::string program;
    std::string description;
    std::vector<argument> args;

    ArgumentParser(std::string&& description);
    void parse_args(const int& argc, char** argv);
    void print_help(const int& error_index);

    void error_parsing_type(const std::string& data_type);
    void error_parsing_inputs_arg();
    void error_parsing_required_arg(const std::string& name);
    void error_parsing_unknown_arg(const std::string& name);
    void error_parsing_invalid_arg(const std::string& name, const std::string& value);
    void error_parsing_out_of_range_arg(const std::string& name, const std::string& value);
    void error_parsing_bool(const std::string& name, const std::string& value);
    void error_parsing_char(const std::string& name, const std::string& value);

    template <typename T>
    void add_argument(std::string&& name, std::string&& short_flag, std::string&& long_flag,
                      std::string&& help, const bool& required, T default_value);

    template <typename T>
    void get(const std::string& name, T& values);
};

template <>
inline void ArgumentParser::add_argument<int>(std::string&& name, std::string&& short_flag, std::string&& long_flag,
                                              std::string&& help, const bool& req, int default_value) {
    std::string type = "int";
    this->args.emplace_back(argument(std::move(name), std::move(short_flag), std::move(long_flag), std::move(help),
                                     req, std::to_string(default_value), std::move(type)));
}

template <>
inline void ArgumentParser::add_argument<long int>(std::string&& name, std::string&& short_flag, std::string&& long_flag,
                                                   std::string&& help, const bool& req, long int default_value) {
    std::string type = "long int";
    this->args.emplace_back(argument(std::move(name), std::move(short_flag), std::move(long_flag), std::move(help),
                                     req, std::to_string(default_value), std::move(type)));
}

template <>
inline void ArgumentParser::add_argument<float>(std::string&& name, std::string&& short_flag, std::string&& long_flag,
                                                std::string&& help, const bool& req, float default_value) {
    std::string type = "float";
    this->args.emplace_back(argument(std::move(name), std::move(short_flag), std::move(long_flag), std::move(help),
                                     req, std::to_string(default_value), std::move(type)));
}

template <>
inline void ArgumentParser::add_argument<double>(std::string&& name, std::string&& short_flag, std::string&& long_flag,
                                                 std::string&& help, const bool& req, double default_value) {
    std::string type = "double";
    this->args.emplace_back(argument(std::move(name), std::move(short_flag), std::move(long_flag), std::move(help),
                                     req, std::to_string(default_value), std::move(type)));
}

template <>
inline void ArgumentParser::add_argument<char>(std::string&& name, std::string&& short_flag, std::string&& long_flag,
                                               std::string&& help, const bool& req, char default_value) {
    std::string type = "char";
    this->args.emplace_back(argument(std::move(name), std::move(short_flag), std::move(long_flag), std::move(help),
                                     req, std::to_string(default_value), std::move(type)));
}

template <>
inline void ArgumentParser::add_argument<bool>(std::string&& name, std::string&& short_flag, std::string&& long_flag,
                                               std::string&& help, const bool& req, bool default_value) {
    std::string type = "bool";
    this->args.emplace_back(argument(std::move(name), std::move(short_flag), std::move(long_flag), std::move(help),
                                     req, std::to_string(default_value), std::move(type)));
}

template <>
inline void ArgumentParser::add_argument<std::string>(std::string&& name, std::string&& short_flag, std::string&& long_flag,
                                                      std::string&& help, const bool& req, std::string default_value) {
    std::string type = "std::string";
    this->args.emplace_back(argument(std::move(name), std::move(short_flag), std::move(long_flag), std::move(help),
                                     req, std::move(default_value), std::move(type)));
}

template <>
inline void ArgumentParser::get<int>(const std::string& name, int& values) {
    for (const auto& arg : args)
        if (arg.name == name) {
            try {
                values = static_cast<int>(std::stod(arg.values[0]));
                return;
            } catch (std::invalid_argument&) {
                this->error_parsing_invalid_arg(name, arg.values[0]);
            } catch (std::out_of_range&) {
                this->error_parsing_out_of_range_arg(name, arg.values[0]);
            }
        }
    this->error_parsing_unknown_arg(name);
}

template <>
inline void ArgumentParser::get<long>(const std::string& name, long& values) {
    for (const auto& arg : args)
        if (arg.name == name) {
            try {
                values = static_cast<long>(std::stod(arg.values[0]));
                return;
            } catch (std::invalid_argument&) {
                this->error_parsing_invalid_arg(name, arg.values[0]);
            } catch (std::out_of_range&) {
                this->error_parsing_out_of_range_arg(name, arg.values[0]);
            }
        }
    this->error_parsing_unknown_arg(name);
}

template <>
inline void ArgumentParser::get<double>(const std::string& name, double& values) {
    for (const auto& arg : args)
        if (arg.name == name) {
            try {
                values = std::stod(arg.values[0]);
                return;
            } catch (std::invalid_argument&) {
                this->error_parsing_invalid_arg(name, arg.values[0]);
            } catch (std::out_of_range&) {
                this->error_parsing_out_of_range_arg(name, arg.values[0]);
            }
        }
    this->error_parsing_unknown_arg(name);
}

template <>
inline void ArgumentParser::get<std::string>(const std::string& name, std::string& values) {
    for (const auto& arg : args)
        if (arg.name == name) {
            values = arg.values[0];
            return;
        }
    this->error_parsing_unknown_arg(name);
}
