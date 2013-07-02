/**

Portions adapted from fec.c.

This work is derived from the "fec" software by Luigi Rizzo, et al., the
copyright notice and licence terms of which are included below for reference.
fec.c -- forward error correction based on Vandermonde matrices 980624 (C)
1997-98 Luigi Rizzo (luigi@iet.unipi.it)

Portions derived from code by Phil Karn (karn@ka9q.ampr.org),
Robert Morelos-Zaragoza (robert@spectra.eng.hawaii.edu) and Hari
Thirumoorthy (harit@spectra.eng.hawaii.edu), Aug 1995

Modifications by Dan Rubenstein
Modifications (C) 1998 Dan Rubenstein (drubenst@cs.umass.edu)


Portions adapted from filefec.py, part of zfec.

Copyright (C) 2007-2010 Allmydata, Inc.
Author: Zooko Wilcox-O'Hearn


Copyright (c) 2013 Richard Kiss

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*/

(function () {

    function log_ceil(n, b) {
        /*
        The smallest integer k such that b^k >= n.
        */
        var p = 1;
        var k = 0;
        while (p < n) {
            p *= b;
            k += 1;
        }
        return k
    }

    function _build_header(m, k, pad, sh) {
        /*
        @param m: the total number of shares; 1 <= m <= 256
        @param k: the number of shares required to reconstruct; 1 <= k <= m
        @param pad: the number of bytes of padding added to the file before encoding; 0 <= pad < k
        @param sh: the shnum of this share; 0 <= k < m

        @return: a compressed string encoding m, k, pad, and sh
        */

        var bitsused = 0;
        var val = 0;

        val |= (m - 1);
        bitsused += 8; // the first 8 bits always encode m

        var kbits = log_ceil(m, 2); // num bits needed to store all possible values of k
        val <<= kbits;
        bitsused += kbits;

        val |= (k - 1);

        var padbits = log_ceil(k, 2); // num bits needed to store all possible values of pad
        val <<= padbits;
        bitsused += padbits;

        val |= pad;

        var shnumbits = log_ceil(m, 2); // num bits needed to store all possible values of shnum
        val <<= shnumbits;
        bitsused += shnumbits;

        val |= sh;

        if (bitsused <= 16) {
            val <<= (16-bitsused);
            return new Uint8Array([(val>>8) & 0xff, (val>>0)&0xff]);
        }
        if (bitsused <= 24) {
            val <<= (24-bitsused);
            return new Uint8Array([(val>>16) & 0xff, (val>>8) & 0xff, (val>>0)&0xff]);
        }
        val <<= (32-bitsused);
        return new Uint8Array([(val>>24) & 0xff, (val>>16) & 0xff, (val>>8) & 0xff, (cs>>0)&0xff]);
    }

    function _parse_header(byte_array) {
        /*
        @param inf: an object which I can call read(1) on to get another byte

        @return: tuple of (m, k, pad, sh,); side-effect: the first one to four
            bytes of inf will be read
        */
        // The first 8 bits always encode m.
        function MASK(bits) {
            return (1<<bits)-1;
        }

        var idx = 2;
        var ch = byte_array[0];
        if (ch===0) {
            throw "Share files corrupted";
        }
        var m = ch + 1;

        // The next few bits encode k.
        var kbits = log_ceil(m, 2); // num bits needed to store all possible values of k
        var b2_bits_left = 8-kbits;
        var kbitmask = MASK(kbits) << b2_bits_left;
        ch = byte_array[1];
        var k = ((ch & kbitmask) >> b2_bits_left) + 1;

        var shbits = log_ceil(m, 2); // num bits needed to store all possible values of shnum
        var padbits = log_ceil(k, 2); // num bits needed to store all possible values of pad

        var val = ch & (~kbitmask);

        var needed_padbits = padbits - b2_bits_left;
        if (needed_padbits > 0) {
            ch = byte_array[idx++];
            val <<= 8;
            val |= ch;
            needed_padbits -= 8;
        }
        //assert needed_padbits <= 0
        var extrabits = -needed_padbits;
        var pad = val >> extrabits;
        val &= MASK(extrabits);

        var needed_shbits = shbits - extrabits;
        if (needed_shbits > 0) {
            ch = byte_array[idx++];
            val <<= 8;
            val |= ch;
            needed_shbits -= 8;
        }
        //assert needed_shbits <= 0

        var gotshbits = -needed_shbits;

        var sh = val >> gotshbits;

        return [m, k, pad, sh, idx];
    }

    var gf_exp = new Uint8Array(510);  /* index->poly form conversion table    */
    var gf_log = new Uint8Array(256); /* Poly->index form conversion table    */
    var inverse = new Uint8Array(256); /* inverse of field elem.               */

    /*
     * modnn(x) computes x % GF_SIZE, where GF_SIZE is 2**GF_BITS - 1,
     * without a slow divide.
     */
    function modnn(x) {
        while (x >= 255) {
            x -= 255;
            x = (x >> 8) + (x & 255);
        }
        return x;
    }

    var gf_mul_table = new Uint8Array(256*256);

    /*
     * Generate GF(2**m) from the irreducible polynomial p(X) in p[0]..p[m]
     * Lookup tables:
     *     index->polynomial form		gf_exp[] contains j= \alpha^i;
     *     polynomial form -> index form	gf_log[ j = \alpha^i ] = i
     * \alpha=x is the primitive element of GF(2^m)
     *
     * For efficiency, gf_exp[] has size 2*GF_SIZE, so that a simple
     * multiplication of two numbers can be resolved without calling modnn
     */

    function _init_mul_table() {
      var i, j;
      for (i = 0; i < 256; i++)
          for (j = 0; j < 256; j++)
              gf_mul_table[i*256+j] = gf_exp[modnn (gf_log[i] + gf_log[j])];

      for (j = 0; j < 256; j++)
          gf_mul_table[0*256+j] = gf_mul_table[j*256+0] = 0;
    }

    _init_mul_table();

    function gf_mul(a, b) {
        return gf_mul_table[a*256+b];
    }

    /*
     * Primitive polynomials - see Lin & Costello, Appendix A,
     * and  Lee & Messerschmitt, p. 453.
     */
    var Pp="101110001";

    function generate_gf() {
        var i;
        var mask;

        mask = 1;                     /* x ** 0 = 1 */
        gf_exp[8] = 0;          /* will be updated at the end of the 1st loop */
        /*
         * first, generate the (polynomial representation of) powers of \alpha,
         * which are stored in gf_exp[i] = \alpha ** i .
         * At the same time build gf_log[gf_exp[i]] = i .
         * The first 8 powers are simply bits shifted to the left.
         */
        for (i = 0; i < 8; i++, mask <<= 1) {
            gf_exp[i] = mask;
            gf_log[gf_exp[i]] = i;
            /*
             * If Pp[i] == 1 then \alpha ** i occurs in poly-repr
             * gf_exp[8] = \alpha ** 8
             */
            if (Pp.charAt(i) === '1')
                gf_exp[8] ^= mask;
        }
        /*
         * now gf_exp[8] = \alpha ** 8 is complete, so can also
         * compute its inverse.
         */
        gf_log[gf_exp[8]] = 8;
        /*
         * Poly-repr of \alpha ** (i+1) is given by poly-repr of
         * \alpha ** i shifted left one-bit and accounting for any
         * \alpha ** 8 term that may occur when poly-repr of
         * \alpha ** i is shifted.
         */
        mask = 1 << 7;
        for (i = 9; i < 255; i++) {
            if (gf_exp[i - 1] >= mask)
                gf_exp[i] = gf_exp[8] ^ ((gf_exp[i - 1] ^ mask) << 1);
            else
                gf_exp[i] = gf_exp[i - 1] << 1;
            gf_log[gf_exp[i]] = i;
        }
        /*
         * log(0) is not defined, so use a special value
         */
        gf_log[0] = 255;
        /* set the extended gf_exp values for fast multiply */
        for (i = 0; i < 255; i++)
            gf_exp[i + 255] = gf_exp[i];

        /*
         * again special cases. 0 has no inverse. This used to
         * be initialized to 255, but it should make no difference
         * since noone is supposed to read from here.
         */
        inverse[0] = 0;
        inverse[1] = 1;
        for (i = 2; i <= 255; i++)
            inverse[i] = gf_exp[255 - gf_log[i]];
    }


    /*
     * Various linear algebra operations that i use often.
     */

    /*
     * addmul() computes dst[] = dst[] + c * src[]
     * This is used often, so better optimize it! Currently the loop is
     * unrolled 16 times, a good value for 486 and pentium-class machines.
     * The case c=0 is also optimized, whereas c=1 is not. These
     * calls are unfrequent in my typical apps so I did not bother.
     */

    function addmul(dst, dst_idx, src, src_idx, c, sz) {
        var i;
        if (c !== 0) {
            for (i=0;i<sz;i++) {
                dst[dst_idx + i] ^= gf_mul(src[src_idx + i], c);
            }
        }
    }

    /*
     * computes C = AB where A is n*k, B is k*m, C is n*m
     */

    function _matmul(a, b, c, n, k, m) {
        var row, col, i;

        for (row = 0; row < n; row++) {
            for (col = 0; col < m; col++) {
                var acc = 0;
                for (i = 0; i < k; i++)
                    acc ^= gf_mul (a[row*k + i], b[col + m*i]);
                c[row * m + col] = acc;
            }
        }
    }

    function memcmp(src, src_idx, dst, dst_idx, size) {
        var i;
        for (i=0;i<size;i++) {
            if (src[src_idx+i] !== dst[dst_idx+i]) {
                return 1;
            }
        }
        return 0;
    }

    /*
     * _invert_mat() takes a matrix and produces its inverse
     * k is the size of the matrix.
     * (Gauss-Jordan, adapted from Numerical Recipes in C)
     * Return non-zero if singular.
     */

    function _invert_mat(src, k) {
        var c, p;
        var irow = 0;
        var icol = 0;
        var row, col, i, ix;

        var indxc = new Uint8Array(k);
        var indxr = new Uint8Array(k);
        var ipiv = new Uint8Array(k);
        var id_row = new Uint8Array(k);

        /*
         * ipiv marks elements already used as pivots.
         */
        for (i = 0; i < k; i++)
            ipiv[i] = 0;

        for (col = 0; col < k; col++) {
            var pivot_row;
            var found_piv = 0;
            /*
             * Zeroing column 'col', look for a non-zero element.
             * First try on the diagonal, if it fails, look elsewhere.
             */
            if (ipiv[col] != 1 && src[col * k + col] != 0) {
                irow = col;
                icol = col;
                found_piv = 1;
            }
            for (row = 0; row < k; row++) {
                if (found_piv) {
                    break;
                }
                if (ipiv[row] != 1) {
                    for (ix = 0; ix < k; ix++) {
                        if (ipiv[ix] === 0) {
                            if (src[row * k + ix] != 0) {
                                irow = row;
                                icol = ix;
                                found_piv = 1;
                                break;
                            }
                        } else
                            {} //assert (ipiv[ix] <= 1);
                    }
                }
            }
          //found_piv:
            ipiv[icol] += 1;
            /*
             * swap rows irow and icol, so afterwards the diagonal
             * element will be correct. Rarely done, not worth
             * optimizing.
             */
            if (irow != icol)
                for (ix = 0; ix < k; ix++) {
                    var tmp = src[irow * k + ix];
                    src[irow * k + ix] = src[icol * k + ix];
                    tmp = src[icol * k + ix];
                }
            indxr[col] = irow;
            indxc[col] = icol;
            //pivot_row = &src[icol * k];
            c = src[icol * k + icol];
            //assert (c != 0);
            if (c != 1) {                       /* otherwhise this is a NOP */
                /*
                 * this is done often , but optimizing is not so
                 * fruitful, at least in the obvious ways (unrolling)
                 */
                c = inverse[c];
                src[icol * k + icol] = 1;
                for (ix = 0; ix < k; ix++)
                    src[icol * k + ix] = gf_mul (c, src[icol * k + ix]);
            }
            /*
             * from all rows, remove multiples of the selected row
             * to zero the relevant entry (in fact, the entry is not zero
             * because we know it must be zero).
             * (Here, if we know that the pivot_row is the identity,
             * we can optimize the addmul).
             */
            id_row[icol] = 1;
            if (memcmp (src, icol*k, id_row, 0, k) !== 0) {
                for (p = 0, ix = 0; ix < k; ix++, p += k) {
                    if (ix !== icol) {
                        c = src[icol + p];
                        src[icol + p] = 0;
                        addmul (src, p, src, icol * k, c, k);
                    }
                }
            }
            id_row[icol] = 0;
        }                           /* done all columns */
        for (col = k; col > 0; col--)
            if (indxr[col-1] != indxc[col-1])
                for (row = 0; row < k; row++) {
                    var tmp = src[row * k + indxr[col-1]];
                    src[row * k + indxr[col-1]] = src[row * k + indxc[col-1]];
                    src[row * k + indxc[col-1]] = tmp;
                }
    }

    /*
     * fast code for inverting a vandermonde matrix.
     *
     * NOTE: It assumes that the matrix is not singular and _IS_ a vandermonde
     * matrix. Only uses the second column of the matrix, containing the p_i's.
     *
     * Algorithm borrowed from "Numerical recipes in C" -- sec.2.8, but largely
     * revised for my purposes.
     * p = coefficients of the matrix (p_i)
     * q = values of the polynomial (known)
     */
    function _invert_vdm (src, k) {
        var i, j, row, col;
        var b, c, p;
        var t, xx;

        if (k === 1)                   /* degenerate case, matrix must be p^0 = 1 */
            return;
        /*
         * c holds the coefficient of P(x) = Prod (x - p_i), i=0..k-1
         * b holds the coefficient for the matrix inversion
         */
        c = new Uint8Array(k);
        b = new Uint8Array(k);

        p = new Uint8Array(k);

        for (j = 1, i = 0; i < k; i++, j += k) {
            c[i] = 0;
            p[i] = src[j];            /* p[i] */
        }
        /*
         * construct coeffs. recursively. We know c[k] = 1 (implicit)
         * and start P_0 = x - p_0, then at each stage multiply by
         * x - p_i generating P_i = x P_{i-1} - p_i P_{i-1}
         * After k steps we are done.
         */
        c[k - 1] = p[0];              /* really -p(0), but x = -x in GF(2^m) */
        for (i = 1; i < k; i++) {
            var p_i = p[i];            /* see above comment */
            for (j = k - 1 - (i - 1); j < k - 1; j++)
                c[j] ^= gf_mul (p_i, c[j + 1]);
            c[k - 1] ^= p_i;
        }

        for (row = 0; row < k; row++) {
            /*
             * synthetic division etc.
             */
            xx = p[row];
            t = 1;
            b[k - 1] = 1;             /* this is in fact c[k] */
            for (i = k - 1; i > 0; i--) {
                b[i-1] = c[i] ^ gf_mul (xx, b[i]);
                t = gf_mul (xx, t) ^ b[i-1];
            }
            for (col = 0; col < k; col++)
                src[col * k + row] = gf_mul (inverse[t], b[col]);
        }
    }

    function init_fec() {
        generate_gf();
        _init_mul_table();
    }

    init_fec();

    function encode(src, block_number_array) {
        var i, j;
        var k;
        var sz = Math.ceil(src.length / this.k);
        var pad = sz*this.k - src.length;
        var fecs = [];
        var header_length;
        if (block_number_array === undefined) {
            block_number_array = [];
            for (i=0;i<this.n;i++) {
                block_number_array.push(i);
            }
        }
        for (i=0;i<block_number_array.length;i++) {
            var header = _build_header(this.n, this.k, pad, i);
            header_length = header.length;
            fecs.push(new Uint8Array(sz+header_length));
            for (j=0;j<header_length;j++) {
                fecs[i][j] = header[j];
            }
        }

        for (k = 0; k < sz; k += this.stride) {
            var stride = ((sz-k) < this.stride)?(sz-k):this.stride;
            for (i=0; i<block_number_array.length; i++) {
                var fecnum = block_number_array[i];
                var fec = fecs[i];
                var p = fecnum * this.k;
                for (j = 0; j < this.k; j++) {
                    addmul(fec, header_length+k, src, j*stride+k*this.k, this.enc_matrix[p+j], stride);
                }
            }
        }
        return fecs;
    }

    function decode(input_array) {
        var row;
        var col;
        var indices = [];
        var idx;
        for (idx=0;idx<input_array.length;idx++) {
            var tuple = _parse_header(input_array[idx]);
            var m, k, pad, sh;
            m = tuple[0];
            k = tuple[1];
            pad = tuple[2];
            indices.push(tuple[3]);
            input_array[idx] = input_array[idx].slice(tuple[4]);
        }

        var size_per_row = input_array[0].length;
        var sz = size_per_row * this.k;
        var output = new Uint8Array(sz);

        for (row=0;row<indices.length;row++) {
            if ((indices[row] < this.k) && (indices[row] != row)) {
                var tmp = input_array[row];
                input_array[row] = input_array[indices[row]];
                input_array[indices[row]] = tmp;
                tmp = indices[row];
                indices[row] = indices[tmp];
                indices[tmp] = tmp;
            }
        }

        /**
         * Build decode matrix into some memory space.
         *
         * @param matrix a space allocated for a k by k matrix
         */
        function build_decode_matrix(fec, indices) {
            var k = fec.k;
            var matrix = new Uint8Array(k * k);
            var i, j;
            var p;
            for (i=0, p=0; i < k; i++, p += k) {
                if (indices[i] < k) {
                    for (j=0;j<k;j++) {
                        matrix[p+j] = 0;
                    }
                    matrix[p+indices[i]] = 1;
                } else {
                    for (j=0;j<k;j++) {
                        matrix[p+j] = fec.enc_matrix[indices[i] * k + j];
                    }
                }
            }
            _invert_mat (matrix, k);
            return matrix;
        }

        var decode_matrix = build_decode_matrix(this, indices);


        var k, col_base = 0;
        for (k = 0; k < size_per_row; k += this.stride) {
            var stride = (size_per_row-k < this.stride)?size_per_row-k:this.stride;
            var out_stride = k*this.k;
            for (row=0; row<this.k; row++) {
                if (indices[row] < this.k) {
                    for (col=0; col<stride; col++) {
                        output[row*stride+out_stride+col] = input_array[indices[row]][k+col];
                    }
                } else {
                    for (col=0; col < this.k; col++) {
                        addmul(output, row*stride+out_stride, input_array[col], k, decode_matrix[row * this.k + col], stride);
                    }
                }
            }
        }
        output = output.subarray(0,sz-pad);
        return output;
    }

    var fec = function(k, n) {
        var row, col;
        var p, tmp_m;

        var enc_matrix = new Uint8Array(n * k);
        tmp_m = new Uint8Array (n * k);
        /*
         * fill the matrix with powers of field elements, starting from 0.
         * The first row is special, cannot be computed with exp. table.
         */
        tmp_m[0] = 1;
        for (col = 1; col < k; col++)
            tmp_m[col] = 0;
        for (p = k, row = 0; row < n - 1; row++, p += k)
            for (col = 0; col < k; col++)
                tmp_m[p + col] = gf_exp[modnn (row * col)];

        /*
         * quick code to build systematic matrix: invert the top
         * k*k vandermonde matrix, multiply right the bottom n-k rows
         * by the inverse, and construct the identity matrix at the top.
         */
        _invert_vdm (tmp_m, k);        /* much faster than _invert_mat */
        _matmul(tmp_m.subarray(k * k), tmp_m, enc_matrix.subarray(k*k, n*k), n - k, k, k);
        /*
         * the upper matrix is I so do not bother with a slow multiply
         */
        var i;
        for (i=0;i<k*k;i++) {
            enc_matrix[i] = 0;
        }
        for (p = 0, col = 0; col < k; col++, p += k + 1) {
            enc_matrix[p] = 1;
        }

        return {
            encode : encode,
            decode : decode,
            k : k,
            n : n,
            enc_matrix : enc_matrix,
            stride : 4096
        }
    }

    var my_module = {
        fec : fec
    }

    // check for node module loader
    if (typeof module !== "undefined" && typeof require !== "undefined") {
        module.exports = my_module;
    } else {
        window["fec"] = my_module;
    }
})();
