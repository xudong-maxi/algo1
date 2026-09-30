/**
 * @file    complex_type.h
 * @brief   Host-test stand-in for the project's complex type used by jdps.h.
 *          In the target project, jdps.h includes the project's own header
 *          instead (see JDPS_COMPLEX_HEADER).
 */
#ifndef COMPLEX_TYPE_H
#define COMPLEX_TYPE_H

typedef struct {
    float r;    /**< real part      */
    float i;    /**< imaginary part */
} complex;

#endif /* COMPLEX_TYPE_H */
