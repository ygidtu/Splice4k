package dsu.progressbar

/**
 * @author Zhangyiming
 * @version 20180902
 * @since 2018.09.02
 *
 * 自定义的进度条
 */

class ProgressBar(private val total: Long? = null, private val message: String = "") {

    private val begin = System.currentTimeMillis()
    private val animationChars = charArrayOf('|', '/', '-', '\\')
    private var current = 0.toLong()
    private var last = this.begin


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

        this.last = tmpTime
        val timeUsage = millisToTime(tmpTime - this.begin)
        if ( total == null ) {
            val speed = this.current.div(timeGap).toInt()
            print("[$timeUsage] ${this.message} - ${animationChars[(this.current % 4).toInt()]} Current: ${this.current} - ${speed}it/s \r")
        } else {
            val progress = (this.current.toDouble() / this.total * 100).toInt()
            val speed = progress.div(timeGap).toInt()
            val eta = millisToTime(((100 - progress) * speed * 1000).toLong())

            var times = progress / 5
            if (times > 20) {
                times = 20
            }

            print("${this.message} - $progress% |${"=".repeat(times)}${" ".repeat(20-times)} | ${speed}it/s [$timeUsage ETA:$eta] \r")
        }
        this.current++
    }

    fun stepTo( step: Long ) {
        this.current = step

        this.step()
    }

    fun clear() {
        this.current = 0.toLong()
    }
}