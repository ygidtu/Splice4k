package dsu.second.index

import dsu.carrier.Exons
import org.apache.log4j.Logger


/**
 * 注释文件index的基础类
 */


abstract class AnnotationIndex(val infile: String) {
    val logger = Logger.getLogger(AnnotationIndex::class.java)

    val data = mutableMapOf<String, MutableList<Exons>>()

    init {
        this.readExon()
    }

    abstract fun getSource(info: List<String>): Map<String, String>

    abstract fun readExon()
}