/**
 * Copyright (c) 2005, www.fontbox.org
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
 * 3. Neither the name of fontbox; nor the names of its
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
 * http://www.fontbox.org
 *
 */
package org.fontbox.ttf;

import java.io.IOException;

import org.fontbox.encoding.MacRomanEncoding;

/**
 * A table in a true type font.
 * 
 * @author Ben Litchfield (ben@benlitchfield.com)
 * @version $Revision: 1.1 $
 */
public class PostScriptTable extends TTFTable
{
    private float formatType;
    private float italicAngle;
    private short underlinePosition;
    private short underlineThickness;
    private long isFixedPitch;
    private long minMemType42;
    private long maxMemType42;
    private long mimMemType1;
    private long maxMemType1;
    private String[] glyphNames = null;
    
    /**
     * A tag that identifies this table type.
     */
    public static final String TAG = "post";
    
    /**
     * This will read the required data from the stream.
     * 
     * @param ttf The font that is being read.
     * @param data The stream to read the data from.
     * @throws IOException If there is an error reading the data.
     */
    public void initData( TrueTypeFont ttf, TTFDataStream data ) throws IOException
    {
        MaximumProfileTable maxp = ttf.getMaximumProfile();
        formatType = data.read32Fixed();
        italicAngle = data.read32Fixed();
        underlinePosition = data.readSignedShort();
        underlineThickness = data.readSignedShort();
        isFixedPitch = data.readUnsignedInt();
        minMemType42 = data.readUnsignedInt();
        maxMemType42 = data.readUnsignedInt();
        mimMemType1 = data.readUnsignedInt();
        maxMemType1 = data.readUnsignedInt();
        MacRomanEncoding encoding = new MacRomanEncoding();
        
        
        if( formatType == 1.0f )
        {
            /*
             * This TrueType font file contains exactly the 258 glyphs in the standard 
             * Macintosh TrueType.
             */
            glyphNames = new String[258];
            for( int i=0; i<glyphNames.length; i++)
            {
                String name = encoding.getName( i );
                if( name != null )
                {
                    glyphNames[i] = name;
                }
            }
        }
        else if( formatType == 2.0f )
        {
            int numGlyphs = data.readUnsignedShort();
            int[] glyphNameIndex = new int[numGlyphs];
            glyphNames = new String[ numGlyphs ];
            int maxIndex = Integer.MIN_VALUE;
            for( int i=0; i<numGlyphs; i++ )
            {
                int index = data.readUnsignedShort();
                glyphNameIndex[i] = index;
                maxIndex = Math.max( maxIndex, index );
            }
            String[] nameArray = null;
            if( maxIndex >= 258 )
            {
                nameArray = new String[ maxIndex-258 +1 ];
                for( int i=0; i<maxIndex-258+1; i++ )
                {
                    int numberOfChars = data.read();
                    nameArray[i]=data.readString( numberOfChars );
                }
            }
            for( int i=0; i<numGlyphs; i++ )
            {
                int index = glyphNameIndex[i];
                if( index < 258 )
                {
                    glyphNames[i] = encoding.getName( index );
                }
                else if( index >= 258 && index <= 32767 )
                {
                    glyphNames[i] = nameArray[index-258];
                }
                else
                {
                    throw new IOException( "Unknown glyph name index:" + index );
                }
            }
        }
        else if( formatType == 2.5f )
        {
            int[] glyphNameIndex = new int[maxp.getNumGlyphs()];
            for( int i=0; i<glyphNameIndex.length; i++)
            {
                int offset = data.readSignedByte();
                glyphNameIndex[i] = i+1+offset;
            }
            glyphNames = new String[glyphNameIndex.length];
            for( int i=0; i<glyphNames.length; i++)
            {
                String name = encoding.getName( glyphNameIndex[i] );
                if( name != null )
                {
                    glyphNames[i] = name;
                }
            }
            
        }
        else if( formatType == 3.0f )
        {
            //no postscript information is provided.
        }
    }
    /**
     * @return Returns the formatType.
     */
    public float getFormatType()
    {
        return formatType;
    }
    /**
     * @param formatTypeValue The formatType to set.
     */
    public void setFormatType(float formatTypeValue)
    {
        this.formatType = formatTypeValue;
    }
    /**
     * @return Returns the isFixedPitch.
     */
    public long getIsFixedPitch()
    {
        return isFixedPitch;
    }
    /**
     * @param isFixedPitchValue The isFixedPitch to set.
     */
    public void setIsFixedPitch(long isFixedPitchValue)
    {
        this.isFixedPitch = isFixedPitchValue;
    }
    /**
     * @return Returns the italicAngle.
     */
    public float getItalicAngle()
    {
        return italicAngle;
    }
    /**
     * @param italicAngleValue The italicAngle to set.
     */
    public void setItalicAngle(float italicAngleValue)
    {
        this.italicAngle = italicAngleValue;
    }
    /**
     * @return Returns the maxMemType1.
     */
    public long getMaxMemType1()
    {
        return maxMemType1;
    }
    /**
     * @param maxMemType1Value The maxMemType1 to set.
     */
    public void setMaxMemType1(long maxMemType1Value)
    {
        this.maxMemType1 = maxMemType1Value;
    }
    /**
     * @return Returns the maxMemType42.
     */
    public long getMaxMemType42()
    {
        return maxMemType42;
    }
    /**
     * @param maxMemType42Value The maxMemType42 to set.
     */
    public void setMaxMemType42(long maxMemType42Value)
    {
        this.maxMemType42 = maxMemType42Value;
    }
    /**
     * @return Returns the mimMemType1.
     */
    public long getMimMemType1()
    {
        return mimMemType1;
    }
    /**
     * @param mimMemType1Value The mimMemType1 to set.
     */
    public void setMimMemType1(long mimMemType1Value)
    {
        this.mimMemType1 = mimMemType1Value;
    }
    /**
     * @return Returns the minMemType42.
     */
    public long getMinMemType42()
    {
        return minMemType42;
    }
    /**
     * @param minMemType42Value The minMemType42 to set.
     */
    public void setMinMemType42(long minMemType42Value)
    {
        this.minMemType42 = minMemType42Value;
    }
    /**
     * @return Returns the underlinePosition.
     */
    public short getUnderlinePosition()
    {
        return underlinePosition;
    }
    /**
     * @param underlinePositionValue The underlinePosition to set.
     */
    public void setUnderlinePosition(short underlinePositionValue)
    {
        this.underlinePosition = underlinePositionValue;
    }
    /**
     * @return Returns the underlineThickness.
     */
    public short getUnderlineThickness()
    {
        return underlineThickness;
    }
    /**
     * @param underlineThicknessValue The underlineThickness to set.
     */
    public void setUnderlineThickness(short underlineThicknessValue)
    {
        this.underlineThickness = underlineThicknessValue;
    }
    /**
     * @return Returns the glyphNames.
     */
    public String[] getGlyphNames()
    {
        return glyphNames;
    }
    /**
     * @param glyphNamesValue The glyphNames to set.
     */
    public void setGlyphNames(String[] glyphNamesValue)
    {
        this.glyphNames = glyphNamesValue;
    }
}
