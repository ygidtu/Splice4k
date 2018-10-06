package com.splice4k.isoforms.base

import com.splice4k.base.Exons
import java.util.Objects

/**
 * @author Zhang Yiming
 * @since 2018.10.05
 * @version 20181005
 * 用于寻找isoforms，构筑加权有向图的边
 */

 class DirectedEdge( private val edge: Pair<Exons, Exons>, var weight: Int = 1 ) {

    fun addWeight( edge: DirectedEdge? ) {
        if ( edge != null ) {
            this.weight += edge.weight
        } else {
            this.weight ++
        }
    }


    /**
     * 获取这条有向边的上游源头
     * @return Exons
     */
    fun from(): Exons {
        return this.edge.first
    }


    /**
     * 获取这条有向边的下游源头
     * @return Exons
     */
    fun to(): Exons {
        return this.edge.second
    }


    override fun hashCode(): Int {
        return Objects.hash( this.edge )
    }

    override fun equals(other: Any?): Boolean {
        if ( other is DirectedEdge ) {
            return this.hashCode() == other.hashCode()
        }
        return false
    }

    override fun toString(): String {
        return "${this.edge.first} -> ${this.edge.second} | ${this.weight}"
    }
}