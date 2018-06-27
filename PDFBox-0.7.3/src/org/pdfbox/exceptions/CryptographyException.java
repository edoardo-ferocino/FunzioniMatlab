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
package org.pdfbox.exceptions;

/**
 * An exception that indicates that something has gone wrong during a
 * cryptography operation.
 *
 * @author <a href="mailto:ben@benlitchfield.com">Ben Litchfield</a>
 * @version $Revision: 1.4 $
 */
public class CryptographyException extends Exception
{
    private Exception embedded;

    /**
     * Constructor.
     *
     * @param msg A msg to go with this exception.
     */
    public CryptographyException( String msg )
    {
        super( msg );
    }

    /**
     * Constructor.
     *
     * @param e The root exception that caused this exception.
     */
    public CryptographyException( Exception e )
    {
        super( e.getMessage() );
        setEmbedded( e );
    }
    /**
     * This will get the exception that caused this exception.
     *
     * @return The embedded exception if one exists.
     */
    public Exception getEmbedded()
    {
        return embedded;
    }
    /**
     * This will set the exception that caused this exception.
     *
     * @param e The sub exception.
     */
    private void setEmbedded( Exception e )
    {
        embedded = e;
    }
}