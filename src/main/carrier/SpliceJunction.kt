package main.carrier

/**
 * @since 2018.06.21
 * @version 0.1
 * @author Zhang Yiming
 *
 * 可变剪接信息的载体
 */

class SpliceJunction(val gene: Genes) {

    data class Event(val name: String, val start: Int, val end: Int) {
        override fun toString(): String {
            return "${this.name}\t${this.start}\t${this.end}"
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
            return "${this.gene.transcriptId}\t${this.gene.start}\t${this.gene.end}\t$event"
        }
        return null
    }
}