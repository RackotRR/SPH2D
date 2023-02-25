#pragma once
#define FULL_STAT_LOGGING

#ifdef FULL_STAT_LOGGING
#include <string>

class PrintLog {
public:
	PrintLog& operator()();
	PrintLog& operator()(const char*);
	PrintLog& operator()(const std::string&);
	PrintLog& operator()(float);
	PrintLog& operator()(long long);
	PrintLog& operator()(int);
	PrintLog& operator()(unsigned);

	static PrintLog& instance() {
		static PrintLog theInstance;
		return theInstance;
	}
	void init(const std::string& experimentName);
	void init_stdout();
private:
	PrintLog() = default;

	void printlog_part(const std::string& line_part);
	void printlog_end_of_line();
};

inline PrintLog& printlog() {
	return PrintLog::instance()();
}
template<typename T>
PrintLog& printlog(T&& value) {
	return PrintLog::instance()(value);
}
inline void init_logger(const std::string& experimentName) {
	PrintLog::instance().init(experimentName);
}
inline void init_logger() {
	PrintLog::instance().init_stdout();
}
#else
struct DontPrintLog {
	DontPrintLog& operator()(...) { return *this; }
};

#define init_logger(...)
#define printlog(...) DontPrintLog()()
#endif // FULL_STAT_LOGGING
