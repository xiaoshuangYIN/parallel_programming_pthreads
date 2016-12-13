#!/bin/bash

./rainfall 10 0.25 4 sample_4x4.in > output_4x4.out
./rainfall 20 0.5 16 sample_16x16.in > output_16x16.out
./rainfall 20 0.5 32 sample_32x32.in > output_32x32.out

diff output_4x4.out sample_4x4.out
diff output_16x16.out sample_16x16.out
diff output_32x32.out sample_32x32.out

rm output_4x4.out
rm output_16x16.out
rm output_32x32.out
