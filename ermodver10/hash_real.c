/*  Written by Shun Sakuraba

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>

uint64_t hash_double(const double *v_, const size_t n)
{
  uint64_t x = 0, buf;
  size_t i;
  for(i = 0; i < n; ++i){
    memcpy(&buf, &v_[i], sizeof(double));
    x = ((x << 7) | (x >> 57)) ^ buf;
  }
  return x;
}

uint64_t hash_float(const float *v_, const size_t n)
{
  uint64_t x = 0;
  uint32_t buf;
  size_t i;
  for(i = 0; i < n; ++i){
    memcpy(&buf, &v_[i], sizeof(float));
    x = ((x << 7) | (x >> 57)) ^ (uint64_t) buf;
  }
  return x;
}
