/*
 * Copyright (c) 1998-2018 University Corporation for Atmospheric Research/Unidata
 * See LICENSE for license information.
 */

package ucar.nc2.iosp.bufr.tables;

import org.jdom2.Element;
import org.jdom2.input.SAXBuilder;
import java.io.InputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Read standard WMO Table A (data categories). */
public class TableA {
  private static org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(TableA.class);
  private static final String TABLEA_FILENAME = "wmo/BUFR_37_0_0_TableA_en.xml";
  private static Map<Integer, String> tableA;

  /*
   * <BUFR_19_1_1_TableA_en>
   * <No>27</No>
   * <CodeFigure>28</CodeFigure>
   * <Meaning_en>Precision orbit (satellite)</Meaning_en>
   * <Status>Operational</Status>
   * </BUFR_19_1_1_TableA_en>
   * 
   * <Exp_BUFRTableA_E>
   * <No>4</No>
   * <CodeFigure>3</CodeFigure>
   * <Meaning_E>Vertical soundings (satellite)</Meaning_E>
   * <Status>Operational</Status>
   * </Exp_BUFRTableA_E>
   */
  private static void init() {
    String filename = BufrTables.RESOURCE_PATH + TABLEA_FILENAME;
    try (InputStream is = CodeFlagTables.class.getResourceAsStream(filename)) {

      HashMap<Integer, String> map = new HashMap<>(100);
      SAXBuilder builder = new SAXBuilder();
      builder.setExpandEntities(false);
      org.jdom2.Document tdoc = builder.build(is);
      org.jdom2.Element root = tdoc.getRootElement();

      List<Element> elems = root.getChildren();
      for (Element elem : elems) {
        String line = elem.getChildText("No");
        String codeS = elem.getChildText("CodeFigure");
        String desc = elem.getChildText("Meaning_en");

        try {
          int code = Integer.parseInt(codeS);
          map.put(code, desc);
        } catch (NumberFormatException e) {
          log.debug("NumberFormatException on line " + line + " in " + codeS);
        }

      }
      tableA = map;

    } catch (Exception e) {
      log.error("Can't read BUFR code table " + filename, e);
    }
  }

  /**
   * data category name, from table A
   *
   * @param cat data category
   * @return category name, or null if not found
   */
  public static String getDataCategory(int cat) {
    if (tableA == null)
      init();
    String result = tableA.get(cat);
    return result != null ? result : "Unknown category=" + cat;
  }

}
