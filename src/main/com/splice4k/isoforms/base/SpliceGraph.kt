package com.splice4k.isoforms.base

import com.splice4k.base.Exons


/**
 * @author Zhang Yiming
 * @since 2018.10.05
 * @version 20181006
 *
 * 用于组装isoforms的剪接图
 */


class SpliceGraph() {
    var nodes: Int = 0
    var edges: Int = 0

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
        this.edges ++
        this.nodes ++
    }


    /**
     * 深度优先搜索
     * 有向图，搜索起来没有特别麻烦
     * 理论上能够获取所有的path
     * @param edge 有向边
     * @param res 收集这条边上可能的exons
     * @return list of exons，路径上涉及到的所有外显子
     */
    private fun depthSearch( edge: DirectedEdge, res: MutableList<Exons> ): MutableList<Exons> {
        if ( edge.to() !in res ) {
            res.add(edge.to())
        }

        if ( this.adj.containsKey(edge.to()) ) {

            for ( i in this.adj[edge.to()]!! ) {
                if ( i.to() in res ) {
                    continue
                }
                this.depthSearch(edge = i, res = res)
            }
        }
        return res
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
                val tmp = mutableListOf(i)
                this.depthSearch( edge = j, res = tmp )

                if ( !results.contains(tmp) ) {
                    results.add(tmp)
                }

            }
        }

        return results
    }

    override fun toString(): String {
        return this.adj.values.map { it.toString() }.joinToString("\n")
    }
}

