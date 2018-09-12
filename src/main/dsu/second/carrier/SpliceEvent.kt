package dsu.second.carrier

import dsu.carrier.Exons
import dsu.carrier.Genes

/**
 * @author Zhangyiming
 * @since 2018.09.05
 * @version 20180905
 * 二代算法记录相应的可变剪接的事件信息的类
 */

class SpliceEvent(val event: String, val subtype: String = "NA") {
    val spliceIn = mutableListOf<Exons>()
    val spliceOut = mutableListOf<Exons>()

    val genes = mutableSetOf<Genes>()


    override fun toString(): String {
        val spliceIn = this.spliceIn.joinToString(separator = ",", prefix = "", postfix = "")
        val spliceOut = this.spliceOut.joinToString(separator = ",", prefix = "", postfix = "")
        var genes = ""
        for ( g in this.genes ) {
            genes += "${g.transcriptId},"
        }
        return "$event\t$subtype\t$spliceIn\t$spliceOut\t$genes"
    }
}
