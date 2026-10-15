/* 主机测试用的 securec 桩（工程中使用真实的 securec.h）。 */
#ifndef SECUREC_H
#define SECUREC_H

#include <string.h>

static inline int memset_s(void* dest, size_t dest_max, int c, size_t count)
{
    if (count > dest_max) {
        return -1;
    }
    memset(dest, c, count);
    return 0;
}

static inline int memcpy_s(void* dest, size_t dest_max, const void* src, size_t count)
{
    if (count > dest_max) {
        return -1;
    }
    memcpy(dest, src, count);
    return 0;
}

#endif  /* SECUREC_H */
