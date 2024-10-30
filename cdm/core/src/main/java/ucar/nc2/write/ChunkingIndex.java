/*
 * Copyright (c) 1998-2020 John Caron and University Corporation for Atmospheric Research/Unidata
 * See LICENSE for license information.
 */
package ucar.nc2.write;

import ucar.ma2.Index;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// contributed by cwardgar@usgs.gov 4/12/2010


/**
 * An index that computes chunk shapes. It is intended to be used to compute the origins and shapes for a series
 * of contiguous writes to a multidimensional array.
 * It writes the first n elements (n < maxChunkElems), then the next, etc.
 * Contributed by cwardgar@usgs.gov 4/12/2010
 * Contributed by mr@mygistar.com 10/30/2024 to support file size > 4G .
 */
public class ChunkingIndex {
  /**
   * Constructor for subclasses only.
   *
   * @param _shape describes an index section: slowest varying comes first (row major)
   */
  public ChunkingIndex(int[] _shape) {
    this.shape = new int[_shape.length]; // optimization over clone
    System.arraycopy(_shape, 0, this.shape, 0, _shape.length);

    rank = shape.length;
    current = new int[rank];
    stride = new int[rank];
    size = computeStrides(shape, stride);
    offset = 0;
    hasvlen = (shape.length > 0 && shape[shape.length - 1] < 0);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  protected int[] shape;
  protected int[] stride;
  protected int rank;

  protected long size; // total number of elements
  protected long offset; // element = offset + stride[0]*current[0] + ...

  protected int[] current; // current element's index, used only for the general case

  protected boolean hasvlen;

  /**
   * Computes the shape of the largest possible <b>contiguous</b> chunk, starting at {@link #getCurrentCounter()}
   * and with {@code numElems <= maxChunkElems}.
   *
   * @param maxChunkElems the maximum number of elements in the chunk shape. The actual element count of the shape
   *        returned is likely to be different, and can be found with {@link Index#computeSize}.
   * @return the shape of the largest possible contiguous chunk.
   */
  public int[] computeChunkShape(long maxChunkElems) {
    int[] chunkShape = new int[rank];

    for (int iDim = 0; iDim < rank; ++iDim) {
      int size = (int) (maxChunkElems / stride[iDim]);
      size = (size == 0) ? 1 : size;
      size = Math.min(size, shape[iDim] - current[iDim]);
      chunkShape[iDim] = size;
    }
    return chunkShape;
  }


  /**
   * Compute total number of elements in the array.
   * Stop at vlen
   *
   * @param shape length of array in each dimension.
   * @return total number of elements in the array.
   */
  public static long computeSize(int[] shape) {
    long product = 1;
    for (int aShape : shape) {
      if (aShape < 0)
        break; // stop at vlen
      product *= aShape;
    }
    return product;
  }

  /**
   * Compute standard strides based on array's shape.
   * Ignore vlen
   *
   * @param shape length of array in each dimension.
   * @param stride put result here
   * @return standard strides based on array's shape.
   */
  private static long computeStrides(int[] shape, int[] stride) {
    long product = 1;
    for (int ii = shape.length - 1; ii >= 0; ii--) {
      int thisDim = shape[ii];
      if (thisDim < 0)
        continue; // ignore vlen
      stride[ii] = (int) product;
      product *= thisDim;
    }
    return product;
  }

  





  /**
   * Get the number of dimensions in the array.
   *
   * @return the number of dimensions in the array.
   */
  public int getRank() {
    return rank;
  }

  /**
   * Get the shape: length of array in each dimension.
   *
   * @return the shape
   */
  public int[] getShape() {
    int[] result = new int[shape.length]; // optimization over clone
    System.arraycopy(shape, 0, result, 0, shape.length);
    return result;
  }

  /**
   * Get the length of the ith dimension.
   *
   * @param index which dimension. must be in [0, getRank())
   * @return the ith dimension length
   */
  public int getShape(int index) {
    return shape[index];
  }

  /**
   * Get the total number of elements in the array.
   *
   * @return the total number of elements in the array.
   */
  public long getSize() {
    return size;
  }

  /**
   * Get the current element's index into the 1D backing array.
   * VLEN stops processing.
   *
   * @return the current element's index into the 1D backing array.
   */
  public long currentElement() {
    long value = offset; // NB: dont have to check each index again
    for (int ii = 0; ii < rank; ii++) { // general rank
      if (shape[ii] < 0)
        break;// vlen
      value += ((long)current[ii]) * stride[ii];
    }
    return value;
  }


  /**
   * Get the current counter.
   *
   * @return copy of the current counter.
   */
  public int[] getCurrentCounter() {
    return current.clone();
  }


  /**
   * Set the current counter from the 1D "current element"
   * currElement = offset + stride[0]*current[0] + ...
   *
   * @param currElement set to this value
   */
  public void setCurrentCounter(long currElement) {
    currElement -= offset;
    for (int ii = 0; ii < rank; ii++) { // general rank
      if (shape[ii] < 0) {
        current[ii] = -1;
        break;
      }
      current[ii] = (int)(currElement / stride[ii]);
      currElement -= ((long)current[ii]) * stride[ii];
    }
    set(current); // transfer to subclass fields
  }

  

  /**
   * Set the current element's index. General-rank case.
   *
   * @param index set current value to these values
   * @return this, so you can use A.get(i.set(i))
   * @throws ArrayIndexOutOfBoundsException if index.length != rank.
   */
  public ChunkingIndex set(int[] index) {
    if (index.length != rank)
      throw new ArrayIndexOutOfBoundsException();
    if (rank == 0)
      return this;
    int prefixrank = (hasvlen ? rank - 1 : rank);
    System.arraycopy(index, 0, current, 0, prefixrank);
    if (hasvlen)
      current[rank] = -1; // ??
    return this;
  }




  public Object clone() {
    ChunkingIndex i;
    try {
      i = (ChunkingIndex) super.clone();
    } catch (CloneNotSupportedException e) {
      return null;
    }
    i.stride = stride.clone();
    i.shape = shape.clone();
    i.current = new int[rank]; // want zeros

    // if (name != null) i.name = name.clone();

    return i;
  }
} 