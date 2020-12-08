/*
 * DoneException.java
 *
 * Created on September 19, 2003, 4:24 PM
 */

/**
 *
 * @author  jknoel
 * @version 
 */
package noel.util.exception;

public class DoneException extends java.lang.Exception {

    /**
     * Creates new <code>OutOfGraphRangeException</code> without detail message.
     */
    public DoneException() {
    }

    /**
     * Constructs an <code>OutOfGraphRangeException</code> with the specified detail message.
     * @param msg the detail message.
     */
    public DoneException(String msg) {
        super(msg);
    }
}