/**
 * Copyright (c) 2003-2006, www.pdfbox.org
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
package org.pdfbox.cos;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.UnsupportedEncodingException;

import org.pdfbox.persistence.util.COSHEXTable;

import org.pdfbox.exceptions.COSVisitorException;

/**
 * This represents a string object in a PDF document.
 *
 * @author <a href="mailto:ben@benlitchfield.com">Ben Litchfield</a>
 * @version $Revision: 1.29 $
 */
public class COSString extends COSBase
{
    /**
     * One of the open string tokens.
     */
    public static final byte[] STRING_OPEN = new byte[]{ 40 }; //"(".getBytes();
    /**
     * One of the close string tokens.
     */
    public static final byte[] STRING_CLOSE = new byte[]{ 41 }; //")".getBytes( "ISO-8859-1" );
    /**
     * One of the open string tokens.
     */
    public static final byte[] HEX_STRING_OPEN = new byte[]{ 60 }; //"<".getBytes( "ISO-8859-1" );
    /**
     * One of the close string tokens.
     */
    public static final byte[] HEX_STRING_CLOSE = new byte[]{ 62 }; //">".getBytes( "ISO-8859-1" );
    /**
     * the escape character in strings.
     */
    public static final byte[] ESCAPE = new byte[]{ 92 }; //"\\".getBytes( "ISO-8859-1" );

    /**
     * CR escape characters.
     */
    public static final byte[] CR_ESCAPE = new byte[]{ 92, 114 }; //"\\r".getBytes( "ISO-8859-1" );
    /**
     * LF escape characters.
     */
    public static final byte[] LF_ESCAPE = new byte[]{ 92, 110 }; //"\\n".getBytes( "ISO-8859-1" );
    /**
     * HT escape characters.
     */
    public static final byte[] HT_ESCAPE = new byte[]{ 92, 116 }; //"\\t".getBytes( "ISO-8859-1" );
    /**
     * BS escape characters.
     */
    public static final byte[] BS_ESCAPE = new byte[]{ 92, 98 }; //"\\b".getBytes( "ISO-8859-1" );
    /**
     * FF escape characters.
     */
    public static final byte[] FF_ESCAPE = new byte[]{ 92, 102 }; //"\\f".getBytes( "ISO-8859-1" );
    
    private ByteArrayOutputStream out = new ByteArrayOutputStream();

    /** 
     * Forces the string to be serialized in literal form but not hexa form. 
     */
    private boolean forceLiteralForm = false;
    
    
    /**
     * Constructor.
     */
    public COSString()
    {
    }

    /**
     * Explicit constructor for ease of manual PDF construction.
     *
     * @param value The string value of the object.
     */
    public COSString( String value )
    {
        try
        {
            boolean unicode16 = false;
            char[] chars = value.toCharArray();
            for( int i=0; i<chars.length; i++ )
            {
                if( chars[i] > 255 )
                {
                    unicode16 = true; 
                }
            }
            if( unicode16 )
            {
                out.write( 0xFE );
                out.write( 0xFF );
                out.write( value.getBytes( "UTF-16BE" ) );
            }
            else
            {
                out.write(value.getBytes("ISO-8859-1"));
            }
        }
        catch (IOException ignore)
        {
            ignore.printStackTrace();
            //should never happen
        }
    }

    /**
     * Explicit constructor for ease of manual PDF construction.
     *
     * @param value The string value of the object.
     */
    public COSString( byte[] value )
    {
        try
        {
            out.write( value );
        }
        catch (IOException ignore)
        {
            ignore.printStackTrace();
            //should never happen
        }
    }
    
    /**
     * Forces the string to be written in literal form instead of hexadecimal form. 
     * 
     * @param v if v is true the string will be written in literal form, otherwise it will
     * be written in hexa if necessary.
     */
    
    public void setForceLiteralForm(boolean v)
    {
        forceLiteralForm = v;
    }
    
    /**
     * This will create a COS string from a string of hex characters.
     * 
     * @param hex A hex string.
     * @return A cos string with the hex characters converted to their actual bytes.
     * @throws IOException If there is an error with the hex string.
     */
    public static COSString createFromHexString( String hex ) throws IOException
    {
        COSString retval = new COSString();
        StringBuffer hexBuffer = new StringBuffer( hex.trim() );
        //if odd number then the last hex digit is assumed to be 0
        if( hexBuffer.length() % 2 == 1 )
        {
            hexBuffer.append( "0" );
        }
        for( int i=0; i<hexBuffer.length();)
        {
            String hexChars = "" + hexBuffer.charAt( i++ ) + hexBuffer.charAt( i++ );
            try
            {
                retval.append( Integer.parseInt( hexChars, 16 ) );
            }
            catch( NumberFormatException e )
            {
                throw new IOException( "Error: Expected hex number, actual='" + hexChars + "'" );
            }
        }
        return retval;
    }
    
    /**
     * This will take this string and create a hex representation of the bytes that make the string.
     * 
     * @return A hex string representing the bytes in this string.
     */
    public String getHexString()
    {
        StringBuffer retval = new StringBuffer( out.size() * 2 );
        byte[] data = getBytes();
        for( int i=0; i<data.length; i++ )
        {
            retval.append( COSHEXTable.HEX_TABLE[ (data[i]+256)%256 ] );
        }
        
        return retval.toString();
    }

    /**
     * This will get the string that this object wraps.
     *
     * @return The wrapped string.
     */
    public String getString()
    {
        String retval;
        String encoding = "ISO-8859-1";
        byte[] data = getBytes();
        int start = 0;
        if( data.length > 2 )
        {
            if( data[0] == (byte)0xFF && data[1] == (byte)0xFE )
            {
                encoding = "UTF-16LE";
                start=2;
            }
            else if( data[0] == (byte)0xFE && data[1] == (byte)0xFF )
            {
                encoding = "UTF-16BE";
                start=2;
            }
        }
        try
        {
            retval = new String( getBytes(), start, data.length-start, encoding );
        }
        catch( UnsupportedEncodingException e )
        {
            //should never happen
            e.printStackTrace();
            retval = new String( getBytes() );
        }
        return retval;
    }

    /**
     * This will append a byte[] to the string.
     *
     * @param data The byte[] to add to this string.
     *
     * @throws IOException If an IO error occurs while writing the byte.
     */
    public void append( byte[] data ) throws IOException
    {
        out.write( data );
    }

    /**
     * This will append a byte to the string.
     *
     * @param in The byte to add to this string.
     *
     * @throws IOException If an IO error occurs while writing the byte.
     */
    public void append( int in ) throws IOException
    {
        out.write( in );
    }

    /**
     * This will reset the internal buffer.
     */
    public void reset()
    {
        out.reset();
    }

    /**
     * This will get the bytes of the string.
     *
     * @return A byte array that represents the string.
     */
    public byte[] getBytes()
    {
        return out.toByteArray();
    }

    /**
     * {@inheritDoc}
     */
    public String toString()
    {
        return "COSString{" + new String( getBytes() ) + "}";
    }
    
    /**
     * This will output this string as a PDF object.
     *  
     * @param output The stream to write to.
     * @throws IOException If there is an error writing to the stream.
     */
    public void writePDF( OutputStream output ) throws IOException
    {
        boolean outsideASCII = false;
        //Lets first check if we need to escape this string.
        byte[] bytes = getBytes();
        for( int i=0; i<bytes.length && !outsideASCII; i++ )
        {
            //if the byte is negative then it is an eight bit byte and is
            //outside the ASCII range.
            outsideASCII = bytes[i] <0;
        }
        if( !outsideASCII || forceLiteralForm )
        {
            output.write(STRING_OPEN);
            for( int i=0; i<bytes.length; i++ )
            {
                int b = (bytes[i]+256)%256;
                switch( b )
                {
                    case '(':
                    case ')':
                    case '\\':
                    {
                        output.write(ESCAPE);
                        output.write(b);
                        break;
                    }
                    case 10: //LF
                    {
                        output.write( LF_ESCAPE );
                        break;
                    }
                    case 13: // CR
                    {
                        output.write( CR_ESCAPE );
                        break;
                    }
                    case '\t':
                    {
                        output.write( HT_ESCAPE );
                        break;
                    }
                    case '\b':
                    {
                        output.write( BS_ESCAPE );
                        break;
                    }
                    case '\f':
                    {
                        output.write( FF_ESCAPE );
                        break;
                    }
                    default:
                    {
                        output.write( b );
                    }
                }
            }
            output.write(STRING_CLOSE);
        }
        else
        {
            output.write(HEX_STRING_OPEN);
            for(int i=0; i<bytes.length; i++ )
            {
                output.write( COSHEXTable.TABLE[ (bytes[i]+256)%256 ] );
            }
            output.write(HEX_STRING_CLOSE);
        }
    }



    /**
     * visitor pattern double dispatch method.
     *
     * @param visitor The object to notify when visiting this object.
     * @return any object, depending on the visitor implementation, or null
     * @throws COSVisitorException If an error occurs while visiting this object.
     */
    public Object accept(ICOSVisitor visitor) throws COSVisitorException
    {
        return visitor.visitFromString( this );
    }

    /**
     * {@inheritDoc}
     */
    public boolean equals(Object obj)
    {
        return (obj instanceof COSString) && java.util.Arrays.equals(((COSString) obj).getBytes(), getBytes());
    }

    /**
     * {@inheritDoc}
     */
    public int hashCode()
    {
        return getBytes().hashCode();
    }
}