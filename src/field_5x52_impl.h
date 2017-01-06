// Copyright (c) 2013 Pieter Wuille
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef _SECP256K1_FIELD_REPR_IMPL_H_
#define _SECP256K1_FIELD_REPR_IMPL_H_

#if defined HAVE_CONFIG_H
#include "libsecp256k1-config.h"
#endif

#include <assert.h>
#include <string.h>
#include "num.h"
#include "field.h"

#if defined(USE_FIELD_5X52_ASM)
#include "field_5x52_asm_impl.h"
#elif defined(USE_FIELD_5X52_INT128)
#include "field_5x52_int128_impl.h"
#else
#error "Please select field_5x52 implementation"
#endif

/** Implements arithmetic modulo FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE FFFFFC2F,
 *  represented as 5 uint64_t's in base 2^52. The values are allowed to contain >52 each. In particular,
 *  each FieldElem has a 'magnitude' associated with it. Internally, a magnitude M means each element
 *  is at most M*(2^53-1), except the most significant one, which is limited to M*(2^49-1). All operations
 *  accept any input with magnitude at most M, and have different rules for propagating magnitude to their
 *  output.
 */

void static secp256k1_fe_inner_start(void) {}
void static secp256k1_fe_inner_stop(void) {}

void static secp256k1_fe_normalize(secp256k1_fe_t *r) {
    uint64_t t0 = r->n[0], t1 = r->n[1], t2 = r->n[2], t3 = r->n[3], t4 = r->n[4];

    // Reduce t4 at the start so there will be at most a single carry from the first pass
    uint64_t x = t4 >> 48; t4 &= 0x0FFFFFFFFFFFFULL;
    uint64_t m;

    // The first pass ensures the magnitude is 1, ...
    t0 += x * 0x1000003D1ULL;
    t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
    t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL; m = t1;
    t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL; m &= t2;
    t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL; m &= t3;

    // ... except for a possible carry at bit 48 of t4 (i.e. bit 256 of the field element)
    assert(t4 >> 49 == 0);

    // At most a single final reduction is needed; check if the value is >= the field characteristic
    x = (t4 >> 48) | ((t4 == 0x0FFFFFFFFFFFFULL) & (m == 0xFFFFFFFFFFFFFULL)
        & (t0 >= 0xFFFFEFFFFFC2FULL));

    // Apply the final reduction (for constant-time behaviour, we do it always)
    t0 += x * 0x1000003D1ULL;
    t1 += (t0 >> 52); t0 &= 0xFFFFFFFFFFFFFULL;
    t2 += (t1 >> 52); t1 &= 0xFFFFFFFFFFFFFULL;
    t3 += (t2 >> 52); t2 &= 0xFFFFFFFFFFFFFULL;
    t4 += (t3 >> 52); t3 &= 0xFFFFFFFFFFFFFULL;

    // If t4 didn't carry to bit 48 already, then it should have after any final reduction
    assert(t4 >> 48 == x);

    // Mask off the possible multiple of 2^256 from the final reduction
    t4 &= 0x0FFFFFFFFFFFFULL;

    r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;

#ifdef VERIFY
    r->magnitude = 1;
    r->normalized = 1;
#endif
}

void static inline secp256k1_fe_set_int(secp256k1_fe_t *r, int a) {
    r->n[0] = a;
    r->n[1] = r->n[2] = r->n[3] = r->n[4] = 0;
#ifdef VERIFY
    r->magnitude = 1;
    r->normalized = 1;
#endif
}

// TODO: not constant time!
int static inline secp256k1_fe_is_zero(const secp256k1_fe_t *a) {
#ifdef VERIFY
    assert(a->normalized);
#endif
    return (a->n[0] == 0 && a->n[1] == 0 && a->n[2] == 0 && a->n[3] == 0 && a->n[4] == 0);
}

int static inline secp256k1_fe_is_odd(const secp256k1_fe_t *a) {
#ifdef VERIFY
    assert(a->normalized);
#endif
    return a->n[0] & 1;
}

// TODO: not constant time!
int static inline secp256k1_fe_equal(const secp256k1_fe_t *a, const secp256k1_fe_t *b) {
#ifdef VERIFY
    assert(a->normalized);
    assert(b->normalized);
#endif
    return (a->n[0] == b->n[0] && a->n[1] == b->n[1] && a->n[2] == b->n[2] && a->n[3] == b->n[3] && a->n[4] == b->n[4]);
}

void static secp256k1_fe_set_b32(secp256k1_fe_t *r, const unsigned char *a) {
    r->n[0] = r->n[1] = r->n[2] = r->n[3] = r->n[4] = 0;
    for (int i=0; i<32; i++) {
        for (int j=0; j<2; j++) {
            int limb = (8*i+4*j)/52;
            int shift = (8*i+4*j)%52;
            r->n[limb] |= (uint64_t)((a[31-i] >> (4*j)) & 0xF) << shift;
        }
    }
#ifdef VERIFY
    r->magnitude = 1;
    r->normalized = 1;
#endif
}

/** Convert a field element to a 32-byte big endian value. Requires the input to be normalized */
void static secp256k1_fe_get_b32(unsigned char *r, const secp256k1_fe_t *a) {
  r[31] = a->n[0] & 0xFF; // i = 0

  r[30] = (a->n[0] >> 8)  & 0xFF; // i = 1
  r[29] = (a->n[0] >> 16) & 0xFF; // i = 2
  r[28] = (a->n[0] >> 24) & 0xFF; // i = 3
  r[27] = (a->n[0] >> 32) & 0xFF; // i = 4
  r[26] = (a->n[0] >> 40) & 0xFF; // i = 5

  r[25] = ((a->n[0] >> 48) & 0xF) // i = 6
        | ((a->n[1] & 0xF) << 4);

  r[24] = (a->n[1] >> 4)  & 0xFF; // i = 7
  r[23] = (a->n[1] >> 12) & 0xFF; // i = 8
  r[22] = (a->n[1] >> 20) & 0xFF; // i = 9
  r[21] = (a->n[1] >> 28) & 0xFF; // i = 10
  r[20] = (a->n[1] >> 36) & 0xFF; // i = 11
  r[19] = (a->n[1] >> 44) & 0xFF; // i = 12

  r[18] = a->n[2] & 0xFF; // i = 13

  r[17] = (a->n[2] >> 8)  & 0xFF; // i = 14
  r[16] = (a->n[2] >> 16) & 0xFF; // i = 15
  r[15] = (a->n[2] >> 24) & 0xFF; // i = 16
  r[14] = (a->n[2] >> 32) & 0xFF; // i = 17
  r[13] = (a->n[2] >> 40) & 0xFF; // i = 18

  r[12] = ((a->n[2] >> 48) & 0xF) // i = 19
        | ((a->n[3] & 0xF) << 4);

  r[11] = (a->n[3] >> 4)  & 0xFF; // i = 20
  r[10] = (a->n[3] >> 12) & 0xFF; // i = 21
  r[9]  = (a->n[3] >> 20) & 0xFF; // i = 22
  r[8]  = (a->n[3] >> 28) & 0xFF; // i = 23
  r[7]  = (a->n[3] >> 36) & 0xFF; // i = 24
  r[6]  = (a->n[3] >> 44) & 0xFF; // i = 25

  r[5] = a->n[4] & 0xFF; // i = 26

  r[4] = (a->n[4] >> 8)  & 0xFF; // i = 27
  r[3] = (a->n[4] >> 16) & 0xFF; // i = 28
  r[2] = (a->n[4] >> 24) & 0xFF; // i = 29
  r[1] = (a->n[4] >> 32) & 0xFF; // i = 30
  r[0] = (a->n[4] >> 40) & 0xFF; // i = 31
}

void static inline secp256k1_fe_negate(secp256k1_fe_t *r, const secp256k1_fe_t *a, int m) {
#ifdef VERIFY
    assert(a->magnitude <= m);
    r->magnitude = m + 1;
    r->normalized = 0;
#endif
    r->n[0] = 0xFFFFEFFFFFC2FULL * (m + 1) - a->n[0];
    r->n[1] = 0xFFFFFFFFFFFFFULL * (m + 1) - a->n[1];
    r->n[2] = 0xFFFFFFFFFFFFFULL * (m + 1) - a->n[2];
    r->n[3] = 0xFFFFFFFFFFFFFULL * (m + 1) - a->n[3];
    r->n[4] = 0x0FFFFFFFFFFFFULL * (m + 1) - a->n[4];
}

void static inline secp256k1_fe_mul_int(secp256k1_fe_t *r, int a) {
#ifdef VERIFY
    r->magnitude *= a;
    r->normalized = 0;
#endif
    r->n[0] *= a;
    r->n[1] *= a;
    r->n[2] *= a;
    r->n[3] *= a;
    r->n[4] *= a;
}

void static inline secp256k1_fe_add(secp256k1_fe_t *r, const secp256k1_fe_t *a) {
#ifdef VERIFY
    r->magnitude += a->magnitude;
    r->normalized = 0;
#endif
    r->n[0] += a->n[0];
    r->n[1] += a->n[1];
    r->n[2] += a->n[2];
    r->n[3] += a->n[3];
    r->n[4] += a->n[4];
}

void static secp256k1_fe_mul(secp256k1_fe_t *r, const secp256k1_fe_t *a, const secp256k1_fe_t *b) {
#ifdef VERIFY
    assert(a->magnitude <= 8);
    assert(b->magnitude <= 8);
    r->magnitude = 1;
    r->normalized = 0;
#endif
    secp256k1_fe_mul_inner(a->n, b->n, r->n);
}

void static secp256k1_fe_sqr(secp256k1_fe_t *r, const secp256k1_fe_t *a) {
#ifdef VERIFY
    assert(a->magnitude <= 8);
    r->magnitude = 1;
    r->normalized = 0;
#endif
    secp256k1_fe_sqr_inner(a->n, r->n);
}

#endif
