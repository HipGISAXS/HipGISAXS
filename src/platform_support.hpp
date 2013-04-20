
#include <sys/sysctl.h>

static inline int sse3_capable() {
	int has_sse3 = 0;
	size_t length = sizeof(has_sse3);
	int error = sysctlbyname("hw.optional.sse3", &has_sse3, &length, NULL, 0);
	if(0 != error) return 0;
	return hasSSE3;
} // sse3_capable()
