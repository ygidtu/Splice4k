package main.extractor

import java.io.File
import java.io.IOException
import java.io.PrintWriter

import org.apache.log4j.Logger

import main.carrier.Genes

open class Extractor(private val silent: Boolean = false) {
    open val logger = Logger.getLogger(Extractor::class.java)
    protected var data = mutableListOf<Genes>().toList()

    /**
     * 保存至文件
     * @outfile 文件路径，提取gff信息保存榆次
     * @silent 减少信息输出
     */
    fun saveTo(outfile: String) {
        if (!this.silent) logger.info("Write information to ${outfile}")

        val outFile = File(outfile).absoluteFile

        var writer = PrintWriter(System.out)

        try{
            if (!outFile.parentFile.exists()) outFile.parentFile.mkdirs()

            writer = PrintWriter(outFile)

            for (i in this.data) {
                writer.println(i.get().joinToString(
                        prefix = "",
                        postfix = "",
                        separator = "\t"
                ))
            }
        } catch (err: IOException) {
            logger.error(err.message)
            for (i in err.stackTrace) {
                logger.error(i)
            }
        } finally {
            writer.close()
        }

    }
}