#include "Logger.h"


#ifdef FULL_STAT_LOGGING
#include <fstream>
#include <iostream>
#include <memory>
#include "CommonIncl.h"

PrintLog& PrintLog::operator()() {
	printlog_end_of_line();
	return *this;
}
PrintLog& PrintLog::operator()(const char* line) {
	printlog_part(line);
	return *this;
}
PrintLog& PrintLog::operator()(const std::string& line) {
	printlog_part(line);
	return *this;
}
PrintLog& PrintLog::operator()(rr_float value) {
	printlog_part(std::to_string(value));
	return *this;
}
PrintLog& PrintLog::operator()(long long value) {
	printlog_part(std::to_string(value));
	return *this;
}
PrintLog& PrintLog::operator()(unsigned value) {
	printlog_part(std::to_string(value));
	return *this;
}
PrintLog& PrintLog::operator()(int value) {
	printlog_part(std::to_string(value));
	return *this;
}
PrintLog& PrintLog::operator()(size_t value) {
	printlog_part(std::to_string(value));
	return *this;
}

static std::unique_ptr<std::ofstream> log_file;
static std::ostream* log_stream;

void PrintLog::init(const std::string& experimentName) {
	std::string filename = experimentName + "\\" + "log.txt";
	log_file = std::make_unique<std::ofstream>(filename);
	log_stream = log_file.get();
}
void PrintLog::init_stdout() {
	log_stream = &std::cout;
}
void PrintLog::printlog_part(const std::string& line_part) {
	*log_stream << line_part;
}
void PrintLog::printlog_end_of_line() {
	*log_stream << std::endl;
}

#endif // FULL_STAT_LOGGING