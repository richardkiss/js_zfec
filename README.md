js_zfec
=======

This is a port of zfec https://pypi.python.org/pypi/zfec, a fast erasure codec
for Python, into JavaScript.


Overview
========

This package performs two operations, encoding and decoding. Encoding takes
some input data and expands its size by producing extra "check blocks", also
called "secondary blocks". Decoding takes some data -- any combination of
blocks of the original data (called "primary blocks") and "secondary blocks",
and produces the original data.

The encoding is parameterized by two integers, K and N. N is the total number
of blocks produced, and K is how many of those blocks are necessary to
reconstruct the original data. N is required to be at least 1 and at most
256, and k is required to be at least 1 and at most N.


Usage
=====

The fec.js file exports the symbol "fec".

Create an instance, passing in K and N.

        // To encode:

        var array_of_uint8
        var f = fec(3, 10);
        var encoded = f.encode(array_of_uint8);

        // encoded now contains 10 arrays of uint8, each just over 1/3rd the size
        // of the original array_of_uint8 array.


        // To decode

        // three arbitrary parts, from above

        var f = fec(3, 10);
        var parts = [encoded[3], encoded[6], encoded[7]];
        var decoded = f.decode(parts);

        // now, decoded will equal the original array_of_uint8 above.

The data produced is compatible with the Python-based command-line tools
zfec and zunfec.
