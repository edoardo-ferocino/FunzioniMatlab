/**
 * Copyright (c) 2003, www.pdfbox.org
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of pdfbox; nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * http://www.pdfbox.org
 *
 */
package org.pdfbox.encoding;

import org.pdfbox.cos.COSBase;
import org.pdfbox.cos.COSName;

/**
 * This is an interface to a text encoder.
 *
 * @author <a href="mailto:ben@benlitchfield.com">Ben Litchfield</a>
 * @version $Revision: 1.10 $
 */
public class MacRomanEncoding extends Encoding
{
    /**
     * Constructor.
     */
    public MacRomanEncoding()
    {
        addCharacterEncoding( 0101, COSName.getPDFName( "A" ) );
        addCharacterEncoding( 0256, COSName.getPDFName( "AE" ) );
        addCharacterEncoding( 0347, COSName.getPDFName( "Aacute" ) );
        addCharacterEncoding( 0345, COSName.getPDFName( "Acircumflex" ) );
        addCharacterEncoding( 0200, COSName.getPDFName( "Adieresis" ) );
        addCharacterEncoding( 0313, COSName.getPDFName( "Agrave" ) );
        addCharacterEncoding( 0201, COSName.getPDFName( "Aring" ) );
        addCharacterEncoding( 0314, COSName.getPDFName( "Atilde" ) );
        addCharacterEncoding( 0102, COSName.getPDFName( "B" ) );
        addCharacterEncoding( 0103, COSName.getPDFName( "C" ) );
        addCharacterEncoding( 0202, COSName.getPDFName( "Ccedilla" ) );
        addCharacterEncoding( 0104, COSName.getPDFName( "D" ) );
        addCharacterEncoding( 0105, COSName.getPDFName( "E" ) );
        addCharacterEncoding( 0203, COSName.getPDFName( "Eacute" ) );
        addCharacterEncoding( 0346, COSName.getPDFName( "Ecircumflex" ) );
        addCharacterEncoding( 0350, COSName.getPDFName( "Edieresis" ) );
        addCharacterEncoding( 0351, COSName.getPDFName( "Egrave" ) );
        addCharacterEncoding( 0106, COSName.getPDFName( "F" ) );
        addCharacterEncoding( 0107, COSName.getPDFName( "G" ) );
        addCharacterEncoding( 0110, COSName.getPDFName( "H" ) );
        addCharacterEncoding( 0111, COSName.getPDFName( "I" ) );
        addCharacterEncoding( 0352, COSName.getPDFName( "Iacute" ) );
        addCharacterEncoding( 0353, COSName.getPDFName( "Icircumflex" ) );
        addCharacterEncoding( 0354, COSName.getPDFName( "Idieresis" ) );
        addCharacterEncoding( 0355, COSName.getPDFName( "Igrave" ) );
        addCharacterEncoding( 0112, COSName.getPDFName( "J" ) );
        addCharacterEncoding( 0113, COSName.getPDFName( "K" ) );
        addCharacterEncoding( 0114, COSName.getPDFName( "L" ) );
        addCharacterEncoding( 0115, COSName.getPDFName( "M" ) );
        addCharacterEncoding( 0116, COSName.getPDFName( "N" ) );
        addCharacterEncoding( 0204, COSName.getPDFName( "Ntilde" ) );
        addCharacterEncoding( 0117, COSName.getPDFName( "O" ) );
        addCharacterEncoding( 0316, COSName.getPDFName( "OE" ) );
        addCharacterEncoding( 0356, COSName.getPDFName( "Oacute" ) );
        addCharacterEncoding( 0357, COSName.getPDFName( "Ocircumflex" ) );
        addCharacterEncoding( 0205, COSName.getPDFName( "Odieresis" ) );
        addCharacterEncoding( 0361, COSName.getPDFName( "Ograve" ) );
        addCharacterEncoding( 0257, COSName.getPDFName( "Oslash" ) );
        addCharacterEncoding( 0315, COSName.getPDFName( "Otilde" ) );
        addCharacterEncoding( 0120, COSName.getPDFName( "P" ) );
        addCharacterEncoding( 0121, COSName.getPDFName( "Q" ) );
        addCharacterEncoding( 0122, COSName.getPDFName( "R" ) );
        addCharacterEncoding( 0123, COSName.getPDFName( "S" ) );
        addCharacterEncoding( 0124, COSName.getPDFName( "T" ) );
        addCharacterEncoding( 0125, COSName.getPDFName( "U" ) );
        addCharacterEncoding( 0362, COSName.getPDFName( "Uacute" ) );
        addCharacterEncoding( 0363, COSName.getPDFName( "Ucircumflex" ) );
        addCharacterEncoding( 0206, COSName.getPDFName( "Udieresis" ) );
        addCharacterEncoding( 0364, COSName.getPDFName( "Ugrave" ) );
        addCharacterEncoding( 0126, COSName.getPDFName( "V" ) );
        addCharacterEncoding( 0127, COSName.getPDFName( "W" ) );
        addCharacterEncoding( 0130, COSName.getPDFName( "X" ) );
        addCharacterEncoding( 0131, COSName.getPDFName( "Y" ) );
        addCharacterEncoding( 0331, COSName.getPDFName( "Ydieresis" ) );
        addCharacterEncoding( 0132, COSName.getPDFName( "Z" ) );
        addCharacterEncoding( 0141, COSName.getPDFName( "a" ) );
        addCharacterEncoding( 0207, COSName.getPDFName( "aacute" ) );
        addCharacterEncoding( 0211, COSName.getPDFName( "acircumflex" ) );
        addCharacterEncoding( 0253, COSName.getPDFName( "acute" ) );
        addCharacterEncoding( 0212, COSName.getPDFName( "adieresis" ) );
        addCharacterEncoding( 0276, COSName.getPDFName( "ae" ) );
        addCharacterEncoding( 0210, COSName.getPDFName( "agrave" ) );
        addCharacterEncoding( 046, COSName.getPDFName( "ampersand" ) );
        addCharacterEncoding( 0214, COSName.getPDFName( "aring" ) );
        addCharacterEncoding( 0136, COSName.getPDFName( "asciicircum" ) );
        addCharacterEncoding( 0176, COSName.getPDFName( "asciitilde" ) );
        addCharacterEncoding( 052, COSName.getPDFName( "asterisk" ) );
        addCharacterEncoding( 0100, COSName.getPDFName( "at" ) );
        addCharacterEncoding( 0213, COSName.getPDFName( "atilde" ) );
        addCharacterEncoding( 0142, COSName.getPDFName( "b" ) );
        addCharacterEncoding( 0134, COSName.getPDFName( "backslash" ) );
        addCharacterEncoding( 0174, COSName.getPDFName( "bar" ) );
        addCharacterEncoding( 0173, COSName.getPDFName( "braceleft" ) );
        addCharacterEncoding( 0175, COSName.getPDFName( "braceright" ) );
        addCharacterEncoding( 0133, COSName.getPDFName( "bracketleft" ) );
        addCharacterEncoding( 0135, COSName.getPDFName( "bracketright" ) );
        addCharacterEncoding( 0371, COSName.getPDFName( "breve" ) );
        addCharacterEncoding( 0245, COSName.getPDFName( "bullet" ) );
        addCharacterEncoding( 0143, COSName.getPDFName( "c" ) );
        addCharacterEncoding( 0377, COSName.getPDFName( "caron" ) );
        addCharacterEncoding( 0215, COSName.getPDFName( "ccedilla" ) );
        addCharacterEncoding( 0374, COSName.getPDFName( "cedilla" ) );
        addCharacterEncoding( 0242, COSName.getPDFName( "cent" ) );
        addCharacterEncoding( 0366, COSName.getPDFName( "circumflex" ) );
        addCharacterEncoding( 072, COSName.getPDFName( "colon" ) );
        addCharacterEncoding( 054, COSName.getPDFName( "comma" ) );
        addCharacterEncoding( 0251, COSName.getPDFName( "copyright" ) );
        addCharacterEncoding( 0333, COSName.getPDFName( "currency1" ) );
        addCharacterEncoding( 0144, COSName.getPDFName( "d" ) );
        addCharacterEncoding( 0240, COSName.getPDFName( "dagger" ) );
        addCharacterEncoding( 0340, COSName.getPDFName( "daggerdbl" ) );
        addCharacterEncoding( 0241, COSName.getPDFName( "degree" ) );
        addCharacterEncoding( 0254, COSName.getPDFName( "dieresis" ) );
        addCharacterEncoding( 0326, COSName.getPDFName( "divide" ) );
        addCharacterEncoding( 044, COSName.getPDFName( "dollar" ) );
        addCharacterEncoding( 0372, COSName.getPDFName( "dotaccent" ) );
        addCharacterEncoding( 0365, COSName.getPDFName( "dotlessi" ) );
        addCharacterEncoding( 0145, COSName.getPDFName( "e" ) );
        addCharacterEncoding( 0216, COSName.getPDFName( "eacute" ) );
        addCharacterEncoding( 0220, COSName.getPDFName( "ecircumflex" ) );
        addCharacterEncoding( 0221, COSName.getPDFName( "edieresis" ) );
        addCharacterEncoding( 0217, COSName.getPDFName( "egrave" ) );
        addCharacterEncoding( 070, COSName.getPDFName( "eight" ) );
        addCharacterEncoding( 0311, COSName.getPDFName( "ellipsis" ) );
        addCharacterEncoding( 0321, COSName.getPDFName( "emdash" ) );
        addCharacterEncoding( 0320, COSName.getPDFName( "endash" ) );
        addCharacterEncoding( 075, COSName.getPDFName( "equal" ) );
        addCharacterEncoding( 041, COSName.getPDFName( "exclam" ) );
        addCharacterEncoding( 0301, COSName.getPDFName( "exclamdown" ) );
        addCharacterEncoding( 0146, COSName.getPDFName( "f" ) );
        addCharacterEncoding( 0336, COSName.getPDFName( "fi" ) );
        addCharacterEncoding( 065, COSName.getPDFName( "five" ) );
        addCharacterEncoding( 0337, COSName.getPDFName( "fl" ) );
        addCharacterEncoding( 0304, COSName.getPDFName( "florin" ) );
        addCharacterEncoding( 064, COSName.getPDFName( "four" ) );
        addCharacterEncoding( 0332, COSName.getPDFName( "fraction" ) );
        addCharacterEncoding( 0147, COSName.getPDFName( "g" ) );
        addCharacterEncoding( 0247, COSName.getPDFName( "germandbls" ) );
        addCharacterEncoding( 0140, COSName.getPDFName( "grave" ) );
        addCharacterEncoding( 076, COSName.getPDFName( "greater" ) );
        addCharacterEncoding( 0307, COSName.getPDFName( "guillemotleft" ) );
        addCharacterEncoding( 0310, COSName.getPDFName( "guillemotright" ) );
        addCharacterEncoding( 0334, COSName.getPDFName( "guilsinglleft" ) );
        addCharacterEncoding( 0335, COSName.getPDFName( "guilsinglright" ) );
        addCharacterEncoding( 0150, COSName.getPDFName( "h" ) );
        addCharacterEncoding( 0375, COSName.getPDFName( "hungarumlaut" ) );
        addCharacterEncoding( 055, COSName.getPDFName( "hyphen" ) );
        addCharacterEncoding( 0151, COSName.getPDFName( "i" ) );
        addCharacterEncoding( 0222, COSName.getPDFName( "iacute" ) );
        addCharacterEncoding( 0224, COSName.getPDFName( "icircumflex" ) );
        addCharacterEncoding( 0225, COSName.getPDFName( "idieresis" ) );
        addCharacterEncoding( 0223, COSName.getPDFName( "igrave" ) );
        addCharacterEncoding( 0152, COSName.getPDFName( "j" ) );
        addCharacterEncoding( 0153, COSName.getPDFName( "k" ) );
        addCharacterEncoding( 0154, COSName.getPDFName( "l" ) );
        addCharacterEncoding( 074, COSName.getPDFName( "less" ) );
        addCharacterEncoding( 0302, COSName.getPDFName( "logicalnot" ) );
        addCharacterEncoding( 0155, COSName.getPDFName( "m" ) );
        addCharacterEncoding( 0370, COSName.getPDFName( "macron" ) );
        addCharacterEncoding( 0265, COSName.getPDFName( "mu" ) );
        addCharacterEncoding( 0156, COSName.getPDFName( "n" ) );
        addCharacterEncoding( 071, COSName.getPDFName( "nine" ) );
        addCharacterEncoding( 0226, COSName.getPDFName( "ntilde" ) );
        addCharacterEncoding( 043, COSName.getPDFName( "numbersign" ) );
        addCharacterEncoding( 0157, COSName.getPDFName( "o" ) );
        addCharacterEncoding( 0227, COSName.getPDFName( "oacute" ) );
        addCharacterEncoding( 0231, COSName.getPDFName( "ocircumflex" ) );
        addCharacterEncoding( 0232, COSName.getPDFName( "odieresis" ) );
        addCharacterEncoding( 0317, COSName.getPDFName( "oe" ) );
        addCharacterEncoding( 0376, COSName.getPDFName( "ogonek" ) );
        addCharacterEncoding( 0230, COSName.getPDFName( "ograve" ) );
        addCharacterEncoding( 061, COSName.getPDFName( "one" ) );
        addCharacterEncoding( 0273, COSName.getPDFName( "ordfeminine" ) );
        addCharacterEncoding( 0274, COSName.getPDFName( "ordmasculine" ) );
        addCharacterEncoding( 0277, COSName.getPDFName( "oslash" ) );
        addCharacterEncoding( 0233, COSName.getPDFName( "otilde" ) );
        addCharacterEncoding( 0160, COSName.getPDFName( "p" ) );
        addCharacterEncoding( 0246, COSName.getPDFName( "paragraph" ) );
        addCharacterEncoding( 050, COSName.getPDFName( "parenleft" ) );
        addCharacterEncoding( 051, COSName.getPDFName( "parenright" ) );
        addCharacterEncoding( 045, COSName.getPDFName( "percent" ) );
        addCharacterEncoding( 056, COSName.getPDFName( "period" ) );
        addCharacterEncoding( 0341, COSName.getPDFName( "periodcentered" ) );
        addCharacterEncoding( 0344, COSName.getPDFName( "perthousand" ) );
        addCharacterEncoding( 053, COSName.getPDFName( "plus" ) );
        addCharacterEncoding( 0261, COSName.getPDFName( "plusminus" ) );
        addCharacterEncoding( 0161, COSName.getPDFName( "q" ) );
        addCharacterEncoding( 077, COSName.getPDFName( "question" ) );
        addCharacterEncoding( 0300, COSName.getPDFName( "questiondown" ) );
        addCharacterEncoding( 042, COSName.getPDFName( "quotedbl" ) );
        addCharacterEncoding( 0343, COSName.getPDFName( "quotedblbase" ) );
        addCharacterEncoding( 0322, COSName.getPDFName( "quotedblleft" ) );
        addCharacterEncoding( 0323, COSName.getPDFName( "quotedblright" ) );
        addCharacterEncoding( 0324, COSName.getPDFName( "quoteleft" ) );
        addCharacterEncoding( 0325, COSName.getPDFName( "quoteright" ) );
        addCharacterEncoding( 0342, COSName.getPDFName( "quotesinglbase" ) );
        addCharacterEncoding( 047, COSName.getPDFName( "quotesingle" ) );
        addCharacterEncoding( 0162, COSName.getPDFName( "r" ) );
        addCharacterEncoding( 0250, COSName.getPDFName( "registered" ) );
        addCharacterEncoding( 0373, COSName.getPDFName( "ring" ) );
        addCharacterEncoding( 0163, COSName.getPDFName( "s" ) );
        addCharacterEncoding( 0244, COSName.getPDFName( "section" ) );
        addCharacterEncoding( 073, COSName.getPDFName( "semicolon" ) );
        addCharacterEncoding( 067, COSName.getPDFName( "seven" ) );
        addCharacterEncoding( 066, COSName.getPDFName( "six" ) );
        addCharacterEncoding( 057, COSName.getPDFName( "slash" ) );
        addCharacterEncoding( 040, COSName.getPDFName( "space" ) );
        addCharacterEncoding( 0243, COSName.getPDFName( "sterling" ) );
        addCharacterEncoding( 0164, COSName.getPDFName( "t" ) );
        addCharacterEncoding( 063, COSName.getPDFName( "three" ) );
        addCharacterEncoding( 0367, COSName.getPDFName( "tilde" ) );
        addCharacterEncoding( 0252, COSName.getPDFName( "trademark" ) );
        addCharacterEncoding( 062, COSName.getPDFName( "two" ) );
        addCharacterEncoding( 0165, COSName.getPDFName( "u" ) );
        addCharacterEncoding( 0234, COSName.getPDFName( "uacute" ) );
        addCharacterEncoding( 0236, COSName.getPDFName( "ucircumflex" ) );
        addCharacterEncoding( 0237, COSName.getPDFName( "udieresis" ) );
        addCharacterEncoding( 0235, COSName.getPDFName( "ugrave" ) );
        addCharacterEncoding( 0137, COSName.getPDFName( "underscore" ) );
        addCharacterEncoding( 0166, COSName.getPDFName( "v" ) );
        addCharacterEncoding( 0167, COSName.getPDFName( "w" ) );
        addCharacterEncoding( 0170, COSName.getPDFName( "x" ) );
        addCharacterEncoding( 0171, COSName.getPDFName( "y" ) );
        addCharacterEncoding( 0330, COSName.getPDFName( "ydieresis" ) );
        addCharacterEncoding( 0264, COSName.getPDFName( "yen" ) );
        addCharacterEncoding( 0172, COSName.getPDFName( "z" ) );
        addCharacterEncoding( 060, COSName.getPDFName( "zero" ) );
    }
    
    /**
     * Convert this standard java object to a COS object.
     *
     * @return The cos object that matches this Java object.
     */
    public COSBase getCOSObject()
    {
        return COSName.MAC_ROMAN_ENCODING;
    }
}