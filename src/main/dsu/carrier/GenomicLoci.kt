package dsu.carrier

import dsu.errors.ChromosomeException
import java.util.*
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

/**
 * @author Zhangyiming
 * @since 2018.09.03
 * @version 20180903
 *
 * 基础类，用来重新调整整个代码结构
 */


open class GenomicLoci: Comparable<GenomicLoci> {
    var chromosome: String = ""

    var start: Int = -1
        set(value) {
            if (value > 0) field = value
            else throw ChromosomeException("$value < 0")
        }

    var end: Int = -1
        set(value) {
            if (value >= this.start && value > 0) field = value
            else throw ChromosomeException("End ($value) < start ${this.start} or $value < 0")
        }

    /**
     * 构造器
     * @param start 起始位点
     * @param end 终止位点
     */
    constructor(start: Int, end: Int) {
        this.start = start
        this.end = end
    }

    /**
     * 构造器
     * @param chromosome 染色体
     * @param start 起始位点
     * @param end 终止位点
     */
    constructor(chromosome: String, start: Int, end: Int): this(start, end) {
        this.chromosome = chromosome
    }

    /**
     * 重载toString
     */
    override fun toString(): String {
        return when(this.chromosome) {
            "" -> "${this.start},${this.end}"
            else -> "${this.chromosome}:${this.start }-${this.end}"
        }
    }

    /**
     * 重载hashCode
     */
    override fun hashCode(): Int {
        return Objects.hash(this.chromosome, this.start, this.end)
    }

    /**
     * 重载equals
     */
    override fun equals(other: Any?): Boolean {
        return other is GenomicLoci && this.hashCode() == other.hashCode()
    }


    /**
     * 重载的compareTo，用来比较位点的上下游
     */
    override fun compareTo(other: GenomicLoci) = compareValuesBy(this, other, {it.chromosome}, {it.start}, {it.end})

    /**
     * 判断该位点是否为另一个位点的上游
     * @param other Genes
     * @return true是上游，false不是
     */
    fun isUpStream(other: GenomicLoci, distanceError: Int = 0) : Boolean {
        if ( this.chromosome != other.chromosome ) {
            return this.chromosome < other.chromosome
        }

        return other.start - this.end >= abs(distanceError)
    }


    /**
     * 判断该位点是否在另一个位点的下游
     * @param other Genes
     * @return true是下游，false不是
     */
    fun isDownStream(other: GenomicLoci, distanceError: Int = 0) : Boolean {
        if ( this.chromosome != other.chromosome ) {
            return this.chromosome > other.chromosome
        }

        return this.start - other.end >= abs(distanceError)
    }


    /**
     * 判断这个位点是否与另外一个位点有重合
     * @param other Genes
     * @return true有重合，false没有重合
     */
    fun isOverlap(other: GenomicLoci): Boolean {
        return this.chromosome == other.chromosome &&
                this.start <= other.end &&
                this.end >= other.start
    }


    /**
     * 返回两个位点的重合bp数
     * @param other Genes
     * @return null 没有重合，int重合的bp数q
     */
    fun getOverlap(other: GenomicLoci): Int? {
        if (!this.isOverlap(other)) {
            return null
        }
        return min(this.end, other.end) - max(this.start, other.start)
    }


    /**
     * 计算两个文件之间的重合程度
     * @param other 第二个位点
     * @param all boolean值； true -> 计算连个位点并集的重合比例；false -> 交集
     * @return double, 两个位点重合的比例，如果<=0，则没有重合
     */
    fun overlapPercent(other: GenomicLoci, all: Boolean = false): Double {

        if (this.chromosome != other.chromosome) {
            return 0.0
        }

        val region = when (all) {
            false ->  (this.end - this.start).toDouble()
            else -> (max(this.end, other.end) - min(this.start, other.start)).toDouble()
        }

        return (min(this.end, other.end) - max(this.start, other.start)) / region * 100
    }

}