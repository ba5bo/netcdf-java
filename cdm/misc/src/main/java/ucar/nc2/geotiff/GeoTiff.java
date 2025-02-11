/*
 * Copyright (c) 1998-2018 John Caron and University Corporation for Atmospheric Research/Unidata
 * See LICENSE for license information.
 */
package ucar.nc2.geotiff;

import java.io.Closeable;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.io.StringWriter;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.nio.ShortBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Formatter;
import java.util.List;
import ucar.nc2.util.CompareNetcdf2;

/**
 * Low level read/write geotiff files.
 *
 * @author John Caron
 * @author Yuan Ho
 */
public class GeoTiff implements Closeable {
  private static final boolean showBytes = false, debugRead = false, debugReadGeoKey = false;
  private static final boolean showHeaderBytes = false;

  private String filename;
  private RandomAccessFile file;
  private FileChannel channel;
  private List<IFDEntry> tags = new ArrayList<>();
  private ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
  private boolean readonly;

  /**
   * Constructor. Does not open or create the file.
   *
   * @param filename pathname of file
   */
  public GeoTiff(String filename) {
    this.filename = filename;
  }

  /**
   * Close the Geotiff file.
   *
   * @throws java.io.IOException on io error
   */
  public void close() throws IOException {
    if (channel != null) {
      if (!readonly) {
        channel.force(true);
        channel.truncate(nextOverflowData);
      }
      channel.close();
    }
    if (file != null)
      file.close();
  }

  /////////////////////////////////////////////////////////////////////////////
  // writing

  private int headerSize = 8;
  private int firstIFD;
  private int lastIFD;
  private int startOverflowData;
  private int nextOverflowData;

  void addTag(IFDEntry ifd) {
    tags.add(ifd);
  }

  List<IFDEntry> getTags() {
    return tags;
  }

  void deleteTag(IFDEntry ifd) {
    tags.remove(ifd);
  }

  void setTransform(double xStart, double yStart, double xInc, double yInc) {
    // tie the raster 0, 0 to xStart, yStart
    addTag(new IFDEntry(Tag.ModelTiepointTag, FieldType.DOUBLE)
        .setValue(new double[] {0.0, 0.0, 0.0, xStart, yStart, 0.0}));

    // define the "affine transformation" : requires grid to be regular (!)
    addTag(new IFDEntry(Tag.ModelPixelScaleTag, FieldType.DOUBLE).setValue(new double[] {xInc, yInc, 0.0}));
  }

  private List<GeoKey> geokeys = new ArrayList<>();

  void addGeoKey(GeoKey geokey) {
    geokeys.add(geokey);
  }

  private void writeGeoKeys() {
    if (geokeys.isEmpty())
      return;

    // count extras
    int extra_chars = 0;
    int extra_ints = 0;
    int extra_doubles = 0;
    for (GeoKey geokey : geokeys) {
      if (geokey.isDouble)
        extra_doubles += geokey.count();
      else if (geokey.isString)
        extra_chars += geokey.valueString().length() + 1;
      else if (geokey.count() > 1)
        extra_ints += geokey.count();
    }

    int n = (geokeys.size() + 1) * 4;
    int[] values = new int[n + extra_ints];
    double[] dvalues = new double[extra_doubles];
    char[] cvalues = new char[extra_chars];
    int icounter = n;
    int dcounter = 0;
    int ccounter = 0;

    values[0] = 1;
    values[1] = 1;
    values[2] = 0;
    values[3] = geokeys.size();
    int count = 4;
    for (GeoKey geokey : geokeys) {
      values[count++] = geokey.tagCode();

      if (geokey.isDouble) {
        values[count++] = Tag.GeoDoubleParamsTag.getCode(); // extra double values here
        values[count++] = geokey.count();
        values[count++] = dcounter;
        for (int k = 0; k < geokey.count(); k++)
          dvalues[dcounter++] = geokey.valueD(k);

      } else if (geokey.isString) {
        String s = geokey.valueString();
        values[count++] = Tag.GeoAsciiParamsTag.getCode(); // extra double values here
        values[count++] = s.length(); // dont include trailing 0 in the count
        values[count++] = ccounter;
        for (int k = 0; k < s.length(); k++)
          cvalues[ccounter++] = s.charAt(k);
        cvalues[ccounter++] = (char) 0;

      } else if (geokey.count() > 1) { // more than one int value
        values[count++] = Tag.GeoKeyDirectoryTag.getCode(); // extra int values here
        values[count++] = geokey.count();
        values[count++] = icounter;
        for (int k = 0; k < geokey.count(); k++)
          values[icounter++] = geokey.value(k);

      } else { // normal case of one int value
        values[count++] = 0;
        values[count++] = 1;
        values[count++] = geokey.value();
      }
    } // loop over geokeys

    addTag(new IFDEntry(Tag.GeoKeyDirectoryTag, FieldType.SHORT).setValue(values));
    if (extra_doubles > 0)
      addTag(new IFDEntry(Tag.GeoDoubleParamsTag, FieldType.DOUBLE).setValue(dvalues));
    if (extra_chars > 0)
      addTag(new IFDEntry(Tag.GeoAsciiParamsTag, FieldType.ASCII).setValue(new String(cvalues)));
  }

  int writeData(byte[] data, int imageNumber) throws IOException {
    if (file == null)
      init();

    if (imageNumber == 1)
      channel.position(headerSize);
    else
      channel.position(nextOverflowData);

    ByteBuffer buffer = ByteBuffer.wrap(data);
    channel.write(buffer);

    if (imageNumber == 1)
      firstIFD = headerSize + data.length;
    else
      firstIFD = data.length + nextOverflowData;

    return nextOverflowData;
  }

  int writeData(short[] data, int imageNumber) throws IOException {
    if (file == null)
      init();

    if (imageNumber == 1)
      channel.position(headerSize);
    else
      channel.position(nextOverflowData);

    // no way around making a copy
    ByteBuffer direct = ByteBuffer.allocateDirect(2 * data.length);
    ShortBuffer buffer = direct.asShortBuffer();
    buffer.put(data);
    // buffer.flip();
    channel.write(direct);

    if (imageNumber == 1)
      firstIFD = headerSize + 2 * data.length;
    else
      firstIFD = 2 * data.length + nextOverflowData;

    return nextOverflowData;
  }

  int writeData(int[] data, int imageNumber) throws IOException {
    if (file == null)
      init();

    if (imageNumber == 1)
      channel.position(headerSize);
    else
      channel.position(nextOverflowData);

    // no way around making a copy
    ByteBuffer direct = ByteBuffer.allocateDirect(4 * data.length);
    IntBuffer buffer = direct.asIntBuffer();
    buffer.put(data);
    // buffer.flip();
    channel.write(direct);

    if (imageNumber == 1)
      firstIFD = headerSize + 4 * data.length;
    else
      firstIFD = 4 * data.length + nextOverflowData;

    return nextOverflowData;
  }

  int writeData(float[] data, int imageNumber) throws IOException {
    if (file == null)
      init();

    if (imageNumber == 1)
      channel.position(headerSize);
    else
      channel.position(nextOverflowData);

    // no way around making a copy
    ByteBuffer direct = ByteBuffer.allocateDirect(4 * data.length);
    FloatBuffer buffer = direct.asFloatBuffer();
    buffer.put(data);
    // ((Buffer) buffer).flip();
    channel.write(direct);

    if (imageNumber == 1)
      firstIFD = headerSize + 4 * data.length;
    else
      firstIFD = 4 * data.length + nextOverflowData;

    return nextOverflowData;
  }

  void writeMetadata(int imageNumber) throws IOException {
    if (file == null)
      init();

    // geokeys all get added at once
    writeGeoKeys();

    // tags gotta be in order
    Collections.sort(tags);
    if (imageNumber == 1) {
      writeHeader(channel);
    } else {
      // now this is not the first image we need to fill the Offset of nextIFD
      channel.position(lastIFD);
      ByteBuffer buffer = ByteBuffer.allocate(4);
      if (debugRead)
        System.out.println("position before writing nextIFD= " + channel.position() + " IFD is " + firstIFD);
      buffer.putInt(firstIFD);
      ((Buffer) buffer).flip();
      channel.write(buffer);
    }
    writeIFD(channel, firstIFD);
  }

  private int writeHeader(FileChannel channel) throws IOException {
    channel.position(0);

    ByteBuffer buffer = ByteBuffer.allocate(8);
    buffer.put((byte) 'M');
    buffer.put((byte) 'M');
    buffer.putShort((short) 42);
    buffer.putInt(firstIFD);

    ((Buffer) buffer).flip();
    channel.write(buffer);

    return firstIFD;
  }

  public void initTags() {
    tags = new ArrayList<>();
    geokeys = new ArrayList<>();
  }

  private void init() throws IOException {
    file = new RandomAccessFile(filename, "rw");
    channel = file.getChannel();
    if (debugRead)
      System.out.println("Opened file to write: '" + filename + "', size=" + channel.size());
    readonly = false;
  }

  private void writeIFD(FileChannel channel, int start) throws IOException {
    channel.position(start);

    ByteBuffer buffer = ByteBuffer.allocate(2);
    int n = tags.size();
    buffer.putShort((short) n);
    ((Buffer) buffer).flip();
    channel.write(buffer);

    start += 2;
    startOverflowData = start + 12 * tags.size() + 4;
    nextOverflowData = startOverflowData;

    for (IFDEntry elem : tags) {
      writeIFDEntry(channel, elem, start);
      start += 12;
    }
    // firstIFD = startOverflowData;
    // position to where the "next IFD" goes
    channel.position(startOverflowData - 4);
    lastIFD = startOverflowData - 4;
    if (debugRead)
      System.out.println("pos before writing nextIFD= " + channel.position());
    buffer = ByteBuffer.allocate(4);
    buffer.putInt(0);
    ((Buffer) buffer).flip();
    channel.write(buffer);
  }

  private void writeIFDEntry(FileChannel channel, IFDEntry ifd, int start) throws IOException {
    channel.position(start);
    ByteBuffer buffer = ByteBuffer.allocate(12);

    buffer.putShort((short) ifd.tag.getCode());
    buffer.putShort((short) ifd.type.code);
    buffer.putInt(ifd.count);

    int size = ifd.count * ifd.type.size;
    if (size <= 4) {
      int done = writeValues(buffer, ifd);
      for (int k = 0; k < 4 - done; k++) // fill out to 4 bytes
        buffer.put((byte) 0);
      ((Buffer) buffer).flip();
      channel.write(buffer);

    } else { // write offset
      buffer.putInt(nextOverflowData);
      ((Buffer) buffer).flip();
      channel.write(buffer);
      // write data
      channel.position(nextOverflowData);
      ByteBuffer vbuffer = ByteBuffer.allocate(size);
      writeValues(vbuffer, ifd);
      ((Buffer) vbuffer).flip();
      channel.write(vbuffer);
      nextOverflowData += size;
    }
  }

  static int writeValues(ByteBuffer buffer, IFDEntry ifd) {
    int done = 0;

    if (ifd.type == FieldType.ASCII) {
      return writeSValue(buffer, ifd);

    } else if (ifd.type == FieldType.RATIONAL) {
      for (int i = 0; i < ifd.count * 2; i++)
        done += writeIntValue(buffer, ifd, ifd.value[i]);

    } else if (ifd.type == FieldType.FLOAT) {
      for (int i = 0; i < ifd.count; i++)
        buffer.putFloat((float) ifd.valueD[i]);
      done += ifd.count * 4;

    } else if (ifd.type == FieldType.DOUBLE) {
      for (int i = 0; i < ifd.count; i++)
        buffer.putDouble(ifd.valueD[i]);
      done += ifd.count * 8;

    } else {
      for (int i = 0; i < ifd.count; i++)
        done += writeIntValue(buffer, ifd, ifd.value[i]);
    }

    return done;
  }

  static int writeIntValue(ByteBuffer buffer, IFDEntry ifd, int v) {
    switch (ifd.type.code) {
      case 1:
      case 2:
      case 6:
      case 7:
        // unsigned byte and ascii
        // signed byte and undefined (usually treated as raw binary)
        buffer.put((byte) v);
        return 1;
      case 3:
      case 8:
        // unsigned and signed short
        buffer.putShort((short) v);
        return 2;
      case 4:
      case 5:
      case 9:
      case 10:
        // unsigned and signed rational and 32-bit integer
        buffer.putInt(v);
        return 4;
    }
    return 0;
  }

  static int writeSValue(ByteBuffer buffer, IFDEntry ifd) {
    buffer.put(ifd.valueS.getBytes(StandardCharsets.UTF_8));
    int size = ifd.valueS.length();
    if ((size & 1) != 0)
      size++; // check if odd
    return size;
  }

  /////////////////////////////////////////////////////////////////////////////
  // reading

  /**
   * Read the geotiff file, using the filename passed in the constructor.
   *
   * @throws IOException on read error
   */
  public void read() throws IOException {
    file = new RandomAccessFile(filename, "r");
    channel = file.getChannel();
    if (debugRead)
      System.out.println("Opened file to read:'" + filename + "', size=" + channel.size());
    readonly = true;

    int nextOffset = readHeader(channel);
    while (nextOffset > 0) {
      nextOffset = readIFD(channel, nextOffset);
      parseGeoInfo();
    }

    // parseGeoInfo();
  }

  IFDEntry findTag(Tag tag) {
    if (tag == null)
      return null;
    for (IFDEntry ifd : tags) {
      if (ifd.tag == tag)
        return ifd;
    }
    return null;
  }

  private int readHeader(FileChannel channel) throws IOException {
    channel.position(0);

    ByteBuffer buffer = ByteBuffer.allocate(8);
    int n = channel.read(buffer);
    assert n == 8;
    ((Buffer) buffer).flip();
    if (showHeaderBytes) {
      printBytes(System.out, "header", buffer, 4);
      buffer.rewind();
    }

    byte b = buffer.get();
    if (b == 73)
      byteOrder = ByteOrder.LITTLE_ENDIAN;
    buffer.order(byteOrder);
    buffer.position(4);
    int firstIFD = buffer.getInt();
    if (debugRead)
      System.out.println(" firstIFD == " + firstIFD);

    return firstIFD;
  }

  private int readIFD(FileChannel channel, int start) throws IOException {
    channel.position(start);

    ByteBuffer buffer = ByteBuffer.allocate(2);
    buffer.order(byteOrder);

    int n = channel.read(buffer);
    assert n == 2;
    ((Buffer) buffer).flip();
    if (showBytes) {
      printBytes(System.out, "IFD", buffer, 2);
      buffer.rewind();
    }
    short nentries = buffer.getShort();
    if (debugRead)
      System.out.println(" nentries = " + nentries);

    start += 2;
    for (int i = 0; i < nentries; i++) {
      IFDEntry ifd = readIFDEntry(channel, start);
      if (debugRead)
        System.out.println(i + " == " + ifd);

      tags.add(ifd);
      start += 12;
    }

    if (debugRead)
      System.out.println(" looking for nextIFD at pos == " + channel.position() + " start = " + start);
    channel.position(start);
    buffer = ByteBuffer.allocate(4);
    buffer.order(byteOrder);
    assert 4 == channel.read(buffer);
    ((Buffer) buffer).flip();
    int nextIFD = buffer.getInt();
    if (debugRead)
      System.out.println(" nextIFD == " + nextIFD);
    return nextIFD;
  }

  private IFDEntry readIFDEntry(FileChannel channel, int start) throws IOException {
    if (debugRead)
      System.out.println("readIFDEntry starting position to " + start);

    channel.position(start);
    ByteBuffer buffer = ByteBuffer.allocate(12);
    buffer.order(byteOrder);
    int n = channel.read(buffer);
    assert n == 12;
    ((Buffer) buffer).flip();
    if (showBytes)
      printBytes(System.out, "IFDEntry bytes", buffer, 12);

    IFDEntry ifd;
    buffer.position(0);
    int code = readUShortValue(buffer);
    Tag tag = Tag.get(code);
    if (tag == null)
      tag = new Tag(code);
    FieldType type = FieldType.get(readUShortValue(buffer));
    int count = buffer.getInt();

    ifd = new IFDEntry(tag, type, count);

    if (ifd.count * ifd.type.size <= 4) {
      readValues(buffer, ifd);
    } else {
      int offset = buffer.getInt();
      if (debugRead)
        System.out.println("position to " + offset);
      channel.position(offset);
      ByteBuffer vbuffer = ByteBuffer.allocate(ifd.count * ifd.type.size);
      vbuffer.order(byteOrder);
      assert ifd.count * ifd.type.size == channel.read(vbuffer);
      ((Buffer) vbuffer).flip();
      readValues(vbuffer, ifd);
    }

    return ifd;
  }

  static void readValues(ByteBuffer buffer, IFDEntry ifd) {

    if (ifd.type == FieldType.ASCII) {
      ifd.valueS = readSValue(buffer, ifd);

    } else if (ifd.type == FieldType.RATIONAL) {
      ifd.value = new int[ifd.count * 2];
      for (int i = 0; i < ifd.count * 2; i++)
        ifd.value[i] = readIntValue(buffer, ifd);

    } else if (ifd.type == FieldType.FLOAT) {
      ifd.valueD = new double[ifd.count];
      for (int i = 0; i < ifd.count; i++)
        ifd.valueD[i] = (double) buffer.getFloat();

    } else if (ifd.type == FieldType.DOUBLE) {
      ifd.valueD = new double[ifd.count];
      for (int i = 0; i < ifd.count; i++)
        ifd.valueD[i] = buffer.getDouble();

    } else {
      ifd.value = new int[ifd.count];
      for (int i = 0; i < ifd.count; i++)
        ifd.value[i] = readIntValue(buffer, ifd);
    }

  }

  static int readIntValue(ByteBuffer buffer, IFDEntry ifd) {
    switch (ifd.type.code) {
      case 1:
      case 2:
        // unsigned byte and unsigned ascii
        return (int) buffer.get() & 0xff;
      case 3:
        // unsigned short
        return readUShortValue(buffer);
      case 4:
      case 5:
        // unsigned rational and unsigned 32-bit integer
        // Yes, this can lead to truncation. This is a bug
        // in the design of the IFDEntry API
        return (int) (buffer.getInt() & 0xffffffffL);
      case 6:
      case 7:
        // signed byte and "undefined" (usually treated as binary data)
        return (int) buffer.get();
      case 8:
        // signed short
        return (int) buffer.getShort();
      case 9:
      case 10:
        // signed rational and signed 32-bit integer
        return buffer.getInt();
    }
    return 0;
  }

  static int readUShortValue(ByteBuffer buffer) {
    return buffer.getShort() & 0xffff;
  }

  static String readSValue(ByteBuffer buffer, IFDEntry ifd) {
    byte[] dst = new byte[ifd.count];
    buffer.get(dst);
    return new String(dst, StandardCharsets.UTF_8);
  }

  private void printBytes(PrintStream ps, String head, ByteBuffer buffer, int n) {
    ps.print(head + " == ");
    for (int i = 0; i < n; i++) {
      byte b = buffer.get();
      int ub = (b < 0) ? b + 256 : b;
      ps.print(ub + "(");
      ps.write(b);
      ps.print(") ");
    }
    ps.println();
  }

  /////////////////////////////////////////////////////////////////////////////
  // geotiff stuff

  private void parseGeoInfo() {
    IFDEntry keyDir = findTag(Tag.GeoKeyDirectoryTag);

    if (null == keyDir)
      return;

    int nkeys = keyDir.value[3];
    if (debugReadGeoKey)
      System.out.println("parseGeoInfo nkeys = " + nkeys + " keyDir= " + keyDir);
    int pos = 4;

    for (int i = 0; i < nkeys; i++) {
      int id = keyDir.value[pos++];
      int location = keyDir.value[pos++];
      int vcount = keyDir.value[pos++];
      int offset = keyDir.value[pos++];

      GeoKey.Tag tag = GeoKey.Tag.getOrMake(id);

      GeoKey key = null;
      if (location == 0) { // simple case
        key = new GeoKey(id, offset);

      } else { // more than one, or non short value
        IFDEntry data = findTag(Tag.get(location));
        if (data == null) {
          System.out.println("********ERROR parseGeoInfo: cant find Tag code = " + location);
        } else if (data.tag == Tag.GeoDoubleParamsTag) { // double params
          double[] dvalue = new double[vcount];
          System.arraycopy(data.valueD, offset, dvalue, 0, vcount);
          key = new GeoKey(tag, dvalue);

        } else if (data.tag == Tag.GeoKeyDirectoryTag) { // int params
          int[] value = new int[vcount];
          System.arraycopy(data.value, offset, value, 0, vcount);
          key = new GeoKey(tag, value);

        } else if (data.tag == Tag.GeoAsciiParamsTag) { // ascii params
          String value = data.valueS.substring(offset, offset + vcount);
          key = new GeoKey(tag, value);
        }

      }

      if (key != null) {
        keyDir.addGeoKey(key);
        if (debugReadGeoKey)
          System.out.println(" yyy  add geokey=" + key);
      }
    }

  }

  /**
   * Write the geotiff Tag information to out.
   */
  public void showInfo(PrintWriter out) {
    out.println("Geotiff file= " + filename);
    for (int i = 0; i < tags.size(); i++) {
      IFDEntry ifd = tags.get(i);
      out.println(i + " IFDEntry == " + ifd);
    }
  }

  /**
   * Write the geotiff Tag information to a String.
   */
  public String showInfo() {
    StringWriter sw = new StringWriter(5000);
    showInfo(new PrintWriter(sw));
    return sw.toString();
  }

  /** @deprecated do not use */
  @Deprecated
  public void compare(GeoTiff other, Formatter f) {
    CompareNetcdf2.compareLists(tags, other.getTags(), f);
  }

  //////////////////////////////////////////////////////////////////////////

  // testing only
  ByteBuffer testReadData(int offset, int size) throws IOException {
    channel.position(offset);
    ByteBuffer buffer = ByteBuffer.allocate(size);
    buffer.order(byteOrder);
    assert size == channel.read(buffer);
    ((Buffer) buffer).flip();
    return buffer;
  }
}
