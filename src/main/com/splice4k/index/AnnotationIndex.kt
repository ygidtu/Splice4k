package splice4k.index

import org.apache.log4j.Logger
import splice4k.base.Exons
import splice4k.base.Genes


/**
 * 注释文件index的基础类
 */


abstract class AnnotationIndex(
        val infile: String,
        val smrt: Boolean = false
) {
    val logger = Logger.getLogger(AnnotationIndex::class.java)

    val data = mutableMapOf<String, MutableList<Exons>>()

    val transcripts: MutableList<Genes> = mutableListOf()

    init {
        this.readExon()
    }

    abstract fun getSource(info: List<String>): Map<String, String>

    abstract fun readExon()
}