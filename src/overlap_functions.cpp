#include <Rcpp.h>
using namespace Rcpp;

// Rasterizes point data and, for each raster row, finds the
// approximate median X ordinate of the overlap zone.
//
// Input is a matrix of X, Y, flightlineID.
//
// Output is a matrix of Xmid, Ymid, Ncells for each
// raster row. Xmid will be NA for any rows where Ncells
// was zero.
//
// [[Rcpp::export]]
NumericMatrix get_overlap_midpoints(NumericMatrix xyf, double res) {
  int npoints = xyf.nrow();

  double x0 = min(xyf.column(0));
  double x1 = max(xyf.column(0));

  double y0 = min(xyf.column(1));
  double y1 = max(xyf.column(1));

  int nr = 1 + (int)((y1 - y0) / res);
  int nc = 1 + (int)((x1 - x0) / res);

  NumericMatrix rasdata(nr, nc);

  int r, c;
  double f, x;
  for (int i = 0; i < npoints; i++) {
    r = (int) ((xyf(i, 1) - y0) / res);
    c = (int) ((xyf(i, 0) - x0) / res);

    f = xyf(i,2);
    x = rasdata(r,c);

    if (x == 0) {
      // First point for this raster cell so record flight line
      rasdata(r,c) = f;
    }
    else if (x > 0 && x != f) {
      // Already one or more points in this cell from another
      // flightline, so flag as an overlap cell by setting
      // value to negative
      rasdata(r,c) = -1;
    }
  }

  // for each raster row, find approximate median X within overlap
  NumericVector overxs(nc);
  NumericMatrix midxy(nr, 3);

  for (r = 0; r < nr; r++) {
    int k = 0;
    for (c = 0; c < nc; c++) {
      if (rasdata(r, c) < 0) {
        overxs(k) = x0 + res*(c + 0.5);
        k++ ;
      }
    }

    double xmid;
    switch(k) {
      case 0:
        xmid = NA_REAL;
        break;
      case 1:
        xmid = overxs(0);
        break;
      default:
        bool even = k % 2 == 0;
        int i = k / 2;
        xmid = even ? (overxs(i) + overxs(i-1)) / 2.0 : overxs(i);
        break;
    }

    midxy(r, 0) = xmid;
    midxy(r, 1) = y0 + res * (r + 0.5);
    midxy(r, 2) = k;
  }

  return midxy;
}

