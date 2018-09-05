package dsu.third.extractor

import java.io.File
import java.io.IOException
import java.io.PrintWriter

import org.apache.log4j.Logger

import dsu.carrier.Genes

open class Extractor(private val silent: Boolean = false): Iterator<Genes?> {
    private val logger = Logger.getLogger(Extractor::class.java)
    var data = mutableListOf<Genes>().toList()

    var totalLine = 0
    get() {
        if (field != this.data.size) {
            return this.data.size
        }
        return field
    }
    var index = 0


    /**
     * 保存至文件
     * @outfile 文件路径，提取gff信息保存榆次
     * @silent 减少信息输出
     */
    fun saveTo(outfile: String) {
        if (!this.silent) logger.info("Write information to $outfile")

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


    /**
     * 遍历内部data时，判断还有没有后续数据
     * @return Boolean, true表明还有后续数据，false表明没有
     */
    override fun hasNext(): Boolean {
        return this.index < this.totalLine - 1
    }


    /**
     * 返回类内部index的下一个数据
     * @return Genes, 若没有下一条，返回null
     */
    override fun next(): Genes? {
        if ( this.index < this.totalLine ) {
            this.index ++
            return this.data[this.index-1]
        }
        return null
    }

    /**
     * 根据index，获取指定位置的数据
     * @param index 获取数据的位置，负值则从后往前取值
     * @return Genes，如果index不在合法范围内，则返回null
     */
    fun get(index: Int): Genes? {
        val tmpIndex = when(index < 0) {
            true -> this.totalLine + index
            false -> index
        }

        if ( tmpIndex < this.totalLine ) {
            return this.data[tmpIndex]
        }
        return null
    }

    /**
     * 对extractor的数据进行排序
     */
    fun sort() {
        this.data = this.data.sorted()
    }
}