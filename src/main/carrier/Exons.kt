package main.carrier

import kotlin.math.*
/**
 * @author zhangyiming
 * @since 2018.07.04
 * @version 0.1
 * 外显子常用的function 三连
 */

/**
 * 检查第一个外显子是否在第二个的上游
 * @param first 外显子位点 Array(start, end)
 * @param second 外显子位点 Array(start, end)
 * @return true 有。false 没有
 */
fun isUpStream(first: Array<Int>, second: Array<Int>): Boolean {
    return first[1] < second[0]
}


/**
 * 检查第一个外显子是否在第二个的下游
 * @param first 外显子位点 Array(start, end)
 * @param second 外显子位点 Array(start, end)
 * @return true 有。false 没有
 */  
fun isDownStream(first: Array<Int>, second: Array<Int>): Boolean {
    return first[0] > second[1]
}


/**
 * 计算两个文件之间的重合程度
 * @param first 第一个位点的坐标，主要是reads上的exon
 * @param second 第二个位点的坐标，主要是基因的intron
 * @return double, 两个位点重合的比例，如果<=0，则没有重合
 */
fun overlapPercent(first: Array<Int>, second: Array<Int>, all: Boolean = false): Double {
    
    if (all) {
        return (min(first[1], second[1]) - max(first[0], second[0])) /
                    (second[1] - second[0]).toDouble() * 100
    } else {
        return (min(first[1], second[1]) - max(first[0], second[0])) /
                    (max(first[1], second[1]) - min(first[0], second[0])).toDouble() * 100
    }
}