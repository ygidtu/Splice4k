package com.splice4k.progressbar

/**
 * @author Zhangyiming
 * @version 20180926
 * @since 2018.09.02
 *
 * 自定义的进度条
 */

class ProgressBar(private val total: Long? = null, private val message: String = "") {

    private val begin = System.currentTimeMillis()
    private val animationChars = charArrayOf('|', '/', '-', '\\')
    private var current = 0.toLong()
    private var last = this.begin
    private var gap = 0.toLong()


    private fun millisToTime( millis: Long): String {
        val seconds = millis / 1000

        val minutes = seconds.div(60)
        val hours = minutes.div(24)
        return "${hours % 24}:${String.format("%02d", minutes % 60)}:${String.format("%02d", seconds % 60)}"
    }

    fun step() {
        val tmpTime = System.currentTimeMillis()
        val timeGap = when {
            tmpTime - this.last < 0.000001 -> 0.000001
            else -> (tmpTime - this.last).toDouble() / 1000
        }

        if (timeGap > 0.12 ) {
            if (this.last > this.begin) {
                print("\r")
            }

            this.last = tmpTime
            val timeUsage = millisToTime(tmpTime - this.begin)
            val speed = (this.current - this.gap).div(timeGap).toInt()
            this.gap = this.current

            if ( total == null ) {
                print("[$timeUsage] ${this.message} - ${animationChars[(this.current % 4).toInt()]} Current: ${this.current} - ${speed}it/s \r")
            } else {
                val progress = (this.current.toDouble() / this.total * 100).toInt()
                val eta = millisToTime(((100 - progress) * speed * 1000).toLong())

                var times = progress / 5
                if (times > 20) {
                    times = 20
                }

                print("${this.message} - $progress% |${"=".repeat(times)}${" ".repeat(20-times)} | ${speed}it/s [$timeUsage ETA:$eta] \r")
            }
        }

        this.current++

    }

    fun stepTo( step: Long ) {
        this.gap = step - this.current
        this.current = step

        this.step()
    }

    fun close() {
        this.total?.let {
            this.current = this.total
            this.step()
        }
        println("\r")
    }
}