package com.splice4k.isoforms.base

import com.splice4k.base.Exons


/**
 * @author Zhang Yiming
 * @since 2018.10.05
 * @version 20181006
 *
 * 用于组装isoforms的剪接图
 */


class SpliceGraph {

    private val adj = mutableMapOf<Exons, MutableList<DirectedEdge>>()
    private val targets = mutableSetOf<Exons>()



    /**
     * 向图中添加边
     * @param edge 边
     * @param weight 权重，边出现的次数
     */
    fun addEdge( edge: Pair<Exons, Exons>, weight: Int = 1 ) {
        val tmp = DirectedEdge(edge = edge, weight = weight)
        if ( this.adj.containsKey(edge.first) ) {

            val elementIndex = this.adj[edge.first]!!.indexOf(tmp)
            if ( elementIndex >= 1 ) {
                this.adj[edge.first]!![elementIndex].addWeight( tmp )
            } else {
                this.adj[edge.first]!!.add(tmp)
            }

        } else {
            this.adj[edge.first] = mutableListOf(tmp)
        }

        this.targets.add( edge.second )
    }


    /**
     * 深度优先搜索
     * 有向图，搜索起来没有特别麻烦
     * 理论上能够获取所有的path
     * @param edge 有向边
     * @param results 收集最终每条完整的path
     * @param currentPath 目前正在走的path
     * @return list of exons，路径上涉及到的所有外显子
     */
    private fun depthFirstSearch( edge: DirectedEdge, results: MutableList<List<Exons>>, currentPath: MutableList<Exons> ): List<Exons> {
        if ( !currentPath.contains(edge.from()) ) {
            currentPath.add(edge.from())
        }

        if ( !currentPath.contains(edge.to()) ) {
            currentPath.add(edge.to())
        }

        if ( this.adj.containsKey(edge.to()) ) {

            for ( i in this.adj[edge.to()]!! ) {
                // make a copy of currentPath, to make sure, this won't collect different path into this list
                val tmpPath = currentPath.asSequence().toMutableList()
                this.depthFirstSearch(edge = i, results = results, currentPath = tmpPath)

                // make sure this full path is visited
                if ( !this.adj.containsKey(tmpPath.last())) {
                    results.add(tmpPath)
                }

            }
        }

        return currentPath
    }


    /**
     * 找出所有的可能的isoforms exon分布
     * @return List of list of exons
     */
    fun getAllIsoforms(): List<List<Exons>> {
        val results = mutableListOf<List<Exons>>()

        if ( this.adj.isEmpty() ) {
            return results
        }

        // isoforms 必定都是从最开头的那个exon起始的，不太可能从中间开始
        for ( i in this.adj.keys ) {

            if ( i in this.targets ) {
                continue
            }

            for ( j in this.adj[i]!! ) {
                this.depthFirstSearch( edge = j, results = results, currentPath = mutableListOf() )
            }
        }

        return results
    }

    override fun toString(): String {
        return this.adj.values.asSequence().map { it.toString() }.joinToString("\n")
    }
}

