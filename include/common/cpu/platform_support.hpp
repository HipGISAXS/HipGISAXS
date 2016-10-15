/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: platform_support.hpp
 *  Created: Apr 20, 2013
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */


#include <sys/sysctl.h>

static inline int sse3_capable() {
	int has_sse3 = 0;
	size_t length = sizeof(has_sse3);
	int error = sysctlbyname("hw.optional.sse3", &has_sse3, &length, NULL, 0);
	if(0 != error) return 0;
	return hasSSE3;
} // sse3_capable()
