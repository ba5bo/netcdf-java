/*
 * Copyright (c) 1998-2018 University Corporation for Atmospheric Research/Unidata
 * See LICENSE for license information.
 */
package ucar.nc2.grib.grib2;

import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;
import ucar.nc2.Attribute;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFiles;
import ucar.nc2.Variable;
import java.io.IOException;

@RunWith(JUnit4.class)
public class TestMRMS {
  static final String testfile =
      "../grib/src/test/data/MRMS_LowLevelCompositeReflectivity_00.50_20141207-072038.grib2.gz";

  static final String testfile24BitPng =
      "../grib/src/test/data/pngEncoding/24-bit/MRMS_FLASH_HP_MAXUNITSTREAMFLOW_00.00_20210615-190000.grib2";

  @Test
  public void checkVariable() throws IOException {
    try (NetcdfFile nc = NetcdfFiles.open(testfile)) {
      Variable var = nc.findVariable("LowLevelCompositeReflectivity_altitude_above_msl");
      Assert.assertNotNull(var);

      Attribute att = var.findAttribute("missing_value");
      Assert.assertNotNull(att);
      Assert.assertEquals(-99., att.getNumericValue().doubleValue(), 1e-6);

      att = var.findAttribute("_FillValue");
      Assert.assertNotNull(att);
      Assert.assertEquals(-999., att.getNumericValue().doubleValue(), 1e-6);
    }
  }

  @Test
  public void checkVariable24bit() throws IOException {
    try (NetcdfFile nc = NetcdfFiles.open(testfile24BitPng)) {
      Variable var = nc.findVariable("FLASH_HP_MAXUNITSTREAMFLOW_surface");
      Assert.assertNotNull(var);

      Attribute att = var.findAttribute("missing_value");
      Assert.assertNotNull(att);
      Assert.assertEquals(-9999., att.getNumericValue().doubleValue(), 1e-6);

      att = var.findAttribute("_FillValue");
      Assert.assertNotNull(att);
      Assert.assertEquals(-9999., att.getNumericValue().doubleValue(), 1e-6);
    }
  }
}
