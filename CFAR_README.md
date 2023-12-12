# CFAR Readme

## 2d CFAR Implementation

2d CFAR is implemented analogous to the 1d implementation shown in the course.
We select a ring region of guard cells around the CUT. Around the guard cells
we select a ring of training cells. We compute the average of the signal
encountered in the training cells. The training cell ring diameter is 3, and the
guard cell ring has a radius of 1.

To compute the filtered output signal, an offset is added to the computed
average and compared to the signal in the CUT. If bigger, the output is 1.0,
else 0.0. The offset is chosen at 6.0 which correspond to a 7.8 STN level. Care
must be taken to convert the signal via db2pow/pow2db functions when adding
the offset.

Since the computation of the filtered signal needs the input of cells in the
training cells around them, we can not compute the signal in the boudary cells.
Instead, boundary cells are set to 0.0. This is achieved by starting with an
output array initialized to 0.0 and by restricting indices to only CUT cells in
the main loop of the CFAR.
