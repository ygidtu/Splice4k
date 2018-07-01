package main.carrier

/**
 * @since 2018.06.21
 * @version 0.1
 * @author Zhang Yiming
 *
 * 可变剪接信息的载体
 */

class SpliceJunction(val gene: Genes) {

    data class Event(
            val name: String,
            val start: Int,
            val end: Int
    ) {
        /**
         * 重载toString
         * @return 名称，start和end以tab分割
         */
        override fun toString(): String {
            return "${this.name}\t${this.start}\t${this.end}"
        }

        /**
         * 提取其中的位置信息
         * @return 按照start和end之间差距的1/20为范围，分别向上下游扩大这么多
         */
        fun getPosition() : String {
//            val gap = (this.end - this.start) / 20
            val gap = 100
            if (this.start >= gap) {
                return "${this.start - gap}-${this.end + gap}"
            }
            return "0-${this.end + gap}"
        }
    }

    private val events = mutableListOf<Event>()

    var index = 0


    /**
     * 向基因中添加新的可变剪接事件
     * @param name 时间名称
     * @param start 发生起始位点
     * @param end 发生终止位点
     */
    fun addEvent(name: String, start: Int, end: Int) {
        this.events.add(Event(name, start, end))
    }

    /**
     * 获取特定splice event所构成的string
     * @param event 特定剪接事件
     * @return 事件字符串，主要包括事件基因，基因位置，和事件位置以及扩展后的位置（方便画图）
     */
    private fun getEventString(event: Event): String {
        val genePosition = "${this.gene.chrom}:${this.gene.start}-${this.gene.end}${this.gene.strand}"
        val eventString = "${event.name}\t${this.gene.chrom}:${event.start}-${event.end}"
        val eventSpan = "${this.gene.chrom}:${event.getPosition()}"
        return "${this.gene.transcriptId}\t$genePosition\t$eventString\t$eventSpan"
    }

    /**
     * 返回基因上各种可变剪接的信息
     */
    fun next(): String? {
        if (this.index < this.events.size) {
            this.index++
            return getEventString(events[this.index-1])
        }
        return null
    }

    /**
     * 重载
     * @return 所有事件构成的字符串，可以直接写到文件，每行一个事件
     */
    override fun toString(): String {
        val tmp = mutableListOf<String>()

        this.events.forEach { tmp.add(this.getEventString(it)) }

        return tmp.joinToString(prefix = "", separator = "\n", postfix = "")
    }
}