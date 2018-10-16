package com.splice4k.tools

import java.io.File
import java.util.Scanner
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.SAMException
import htsjdk.samtools.SAMFormatException
import java.io.IOException
import java.lang.NumberFormatException

/**
 * 检查文件格式
 * @since 2018.09.28
 * @author Zhang Yiming
 * @version 20180928
 */

class FileValidator() {

    /**
     * 检查是否为BAM/SAM格式的文件
     * @param infile 输入文件
     */
    private fun isBam(infile: File): Boolean {
        try{

            val reader = SamReaderFactory
                    .makeDefault()
                    .open(infile)
                    .iterator()
            for ( i in reader ) {
                break
            }
            reader.close()
            return true
        } catch (e: SAMException) {

        } catch (e: SAMFormatException) {

        } catch (e: java.lang.RuntimeException) {

        }
        return false
    }


    /**
     * 检查是否为sj文件
     * @param infile 输入文件
     */
    private fun isSJ(infile: File): Boolean {
        try{
            val reader = Scanner(infile)
            var res = false
            reader.use {

                while ( reader.hasNext() ) {
                    val line = reader.nextLine()

                    if ( line.startsWith("#") ) {
                        continue
                    }

                    if ( line.matches("^[\\w\\.]+:\\d+-\\d+[+-\\.]?\t\\d+$".toRegex()) ) {
                        res = true
                    }

                    break
                }
            }

            return res
        } catch ( e: IOException ) {

        }
        return false
    }


    /**
     * 检查是否为sj文件
     * @param infile 输入文件
     */
    private fun isStar(infile: File): Boolean {

        try{
            val reader = Scanner(infile.absoluteFile)

            var res = false

            reader.use {
                while ( reader.hasNext() ) {
                    val line = reader.nextLine()

                    if ( line.startsWith("#") ) {
                        continue
                    }

                    // ^([\w\.]+)(\s+\d+){8}(\s+)?
                    if ( line.matches("^([\\w\\.]+)\\s+\\d+\\s+\\d+\\s+[012]\\s+\\d+\\s+[01]\\s+\\d+\\s+\\d+\\s+\\d+$".toRegex()) ) {
                        res = true
                    }

                    break
                }
            }

            return res
        } catch ( e: IOException ) {

        }
        return false
    }


    /**
     * 检查是否为gff文件
     * @param infile 输入文件
     */
    private fun isGff(infile: File): Boolean {
        try{
            val reader = Scanner(infile)

            var res = false

            reader.use {
                while ( reader.hasNext() ) {
                    val line = reader.nextLine()

                    if ( line.startsWith("#") ) {
                        continue
                    }

                    val lines = line.split("\t")

                    if ( lines[8].matches("([\\w-\\.]+=[\\w:\\s-%,\\.]+;)+([\\w-]+=[\\w:\\s-%,\\.]+)?\$".toRegex()) ) {
                        res = true
                    }

                    break
                }
            }

            return res
        } catch ( e: Exception ) {

        }
        return false
    }


    /**
     * 检查是否为gff文件
     * @param infile 输入文件
     */
    private fun isGtf(infile: File): Boolean {
        try{
            val reader = Scanner(infile)

            var res = false

            reader.use {
                while ( reader.hasNext() ) {
                    val line = reader.nextLine()

                    if ( line.startsWith("#") ) {
                        continue
                    }

                    val lines = line.split("\t")

                    if ( lines[8].matches("([\\w-]+ \"[\\w+\\.\\s-%,:]+\"; ?)+".toRegex()) ) {
                        res = true
                    }

                    break
                }
            }

            return res
        } catch ( e: Exception ) {

        }
        return false
    }


    /**
     * 是否为gmap的align输出文件
     */
    private fun isGmap( infile: File ): Boolean {

        return false
    }

    /**
     * 调用的api
     * @param infile 输入文件
     */
    fun check(infile: File): String {
        if ( this.isBam(infile) ) {
            return "bam"
        }

        if ( this.isStar(infile) ) {
            return "star"
        }

        if ( this.isSJ(infile) ) {
            return "sj"
        }

        if ( this.isGff(infile) ) {
            return "gff"
        }

        if ( this.isGtf(infile) ) {
            return "gtf"
        }

        return "gmap"
    }
}
