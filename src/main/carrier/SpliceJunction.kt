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
            val gap = (this.end - this.start) / 20

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
     */
    fun addEvent(name: String, start: Int, end: Int) {
        this.events.add(Event(name, start, end))
    }


    /**
     * 返回基因上各种可变剪接的信息
     */
    fun next(): String? {
        if (this.index < this.events.size) {
            val event = events[this.index]
            this.index++
            val eventString = "${event.name}\t${this.gene.chrom}:${event.getPosition()}"
            return "${this.gene.transcriptId}\t${this.gene.chrom}\t${this.gene.start}\t${this.gene.end}\t$eventString"
        }
        return null
    }
}