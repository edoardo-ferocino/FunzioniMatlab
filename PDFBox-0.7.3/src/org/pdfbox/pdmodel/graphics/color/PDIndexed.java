/**
 * Copyright (c) 2004, www.pdfbox.org
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
package org.pdfbox.pdmodel.graphics.color;

import org.pdfbox.cos.COSArray;
import org.pdfbox.cos.COSBase;
import org.pdfbox.cos.COSInteger;
import org.pdfbox.cos.COSName;
import org.pdfbox.cos.COSNumber;
import org.pdfbox.cos.COSStream;
import org.pdfbox.cos.COSString;

import java.awt.color.ColorSpace;
import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;

import java.io.InputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

/**
 * This class represents an Indexed color space.
 *
 * @author <a href="mailto:ben@benlitchfield.com">Ben Litchfield</a>
 * @version $Revision: 1.4 $
 */
public class PDIndexed extends PDColorSpace
{

    /**
     * The name of this color space.
     */
    public static final String NAME = "Indexed";

    /**
     * The abbreviated name of this color space.
     */
    public static final String ABBREVIATED_NAME = "I";

    private COSArray array;

    /**
     * Constructor, default DeviceRGB, hival 255.
     */
    public PDIndexed()
    {
        array = new COSArray();
        array.add( COSName.getPDFName( NAME ) );
        array.add( COSName.getPDFName( PDDeviceRGB.NAME ) );
        array.add( new COSInteger( 255 ) );
        array.add( org.pdfbox.cos.COSNull.NULL );
    }

    /**
     * Constructor.
     *
     * @param indexedArray The array containing the indexed parameters
     */
    public PDIndexed( COSArray indexedArray )
    {
        array = indexedArray;
    }

    /**
     * This will return the number of color components.  This will return the
     * number of color components in the base color.
     *
     * @return The number of components in this color space.
     *
     * @throws IOException If there is an error getting the number of color components.
     */
    public int getNumberOfComponents() throws IOException
    {
        return getBaseColorSpace().getNumberOfComponents();
    }

    /**
     * This will return the name of the color space.
     *
     * @return The name of the color space.
     */
    public String getName()
    {
        return NAME;
    }

    /**
     * Create a Java colorspace for this colorspace.
     *
     * @return A color space that can be used for Java AWT operations.
     *
     * @throws IOException If there is an error creating the color space.
     */
    public ColorSpace createColorSpace() throws IOException
    {
        throw new IOException( "Not implemented" );
    }
    
    /**
     * Create a Java color model for this colorspace.
     *
     * @param bpc The number of bits per component.
     *
     * @return A color model that can be used for Java AWT operations.
     *
     * @throws IOException If there is an error creating the color model.
     */
    public ColorModel createColorModel( int bpc ) throws IOException
    {
        int size = getHighValue();
        byte[] index = getLookupData();
        //for (int i=0;i<index.length;i++) System.out.print(index[i]+" ");
        
        ColorModel cm = new IndexColorModel(bpc, size+1, index,0,false);
        return cm;
    }

    /**
     * This will get the color space that acts as the index for this color space.
     *
     * @return The base color space.
     *
     * @throws IOException If there is error creating the base color space.
     */
    public PDColorSpace getBaseColorSpace() throws IOException
    {
        PDColorSpace retval = null;
        COSBase base = array.getObject( 1 );
        if( base instanceof COSName )
        {
            retval = PDColorSpaceFactory.createColorSpace( base );
        }
        else
        {
            throw new IOException( "Error:unknown base colorspace" );
        }

        return retval;
    }

    /**
     * This will set the base color space.
     *
     * @param base The base color space to use as the index.
     */
    public void setBaseColorSpace( PDColorSpace base )
    {
        array.set( 1, base.getCOSObject() );
    }

    /**
     * Get the highest value for the lookup.
     *
     * @return The hival entry.
     */
    public int getHighValue()
    {
        return ((COSNumber)array.getObject( 2 )).intValue();
    }

    /**
     * This will set the highest value that is allowed.  This cannot be higher
     * than 255.
     *
     * @param high The highest value for the lookup table.
     */
    public void setHighValue( int high )
    {
        array.set( 2, new COSInteger( high ) );
    }

    /**
     * This will perform a lookup into the color lookup table.
     *
     * @param componentNumber The component number, probably 1,2,3,3.
     * @param lookupIndex The zero-based index into the table, should not exceed the high value.
     *
     * @return The value that was from the lookup table.
     *
     * @throws IOException If there is an error looking up the color.
     */
    public int lookupColor( int componentNumber, int lookupIndex ) throws IOException
    {
        PDColorSpace baseColor = getBaseColorSpace();
        byte[] data = getLookupData();
        int numberOfComponents = baseColor.getNumberOfComponents();
        return (data[componentNumber*numberOfComponents + lookupIndex]+256)%256;
    }

    private byte[] getLookupData() throws IOException
    {
        COSBase lookupTable = array.getObject( 3 );
        byte[] data = null;
        if( lookupTable instanceof COSString )
        {
            data = ((COSString)lookupTable).getBytes();
        }
        else if( lookupTable instanceof COSStream )
        {
            //Data will be small so just load the whole thing into memory for
            //easier processing
            COSStream lookupStream = (COSStream)lookupTable;
            InputStream input = lookupStream.getUnfilteredStream();
            ByteArrayOutputStream output = new ByteArrayOutputStream(1024);
            byte[] buffer = new byte[ 1024 ];
            int amountRead;
            while( (amountRead = input.read(buffer, 0, buffer.length)) != -1 )
            {
                output.write( buffer, 0, amountRead );
            }
            data = output.toByteArray();
        }
        else if( lookupTable == null )
        {
            data = new byte[0];
        }
        else
        {
            throw new IOException( "Error: Unknown type for lookup table " + lookupTable );
        }
        return data;
    }

    /**
     * This will set a color in the color lookup table.
     *
     * @param componentNumber The component number, probably 1,2,3,3.
     * @param lookupIndex The zero-based index into the table, should not exceed the high value.
     * @param color The color that will go into the table.
     *
     * @throws IOException If there is an error looking up the color.
     */
    public void setLookupColor( int componentNumber, int lookupIndex, int color ) throws IOException
    {
        PDColorSpace baseColor = getBaseColorSpace();
        int numberOfComponents = baseColor.getNumberOfComponents();
        byte[] data = getLookupData();
        data[componentNumber*numberOfComponents + lookupIndex] = (byte)color;
        COSString string = new COSString( data );
        array.set( 3, string );
    }
}