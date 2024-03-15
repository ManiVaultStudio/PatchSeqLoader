#include "Taxonomy.h"

#include "json.hpp"

#include <fstream>
#include <iostream>
#include <QDebug>

using json = nlohmann::json;

namespace
{
    std::vector<TaxonomyNode> traverseJsonArray(Taxonomy& taxonomy, json element);
    std::vector<TaxonomyNode> traverseJsonObject(Taxonomy& taxonomy, json element);

    TaxonomyNode parseNodeAttributes(Taxonomy& taxonomy, json element)
    {
        Q_ASSERT(element.size() == 1);

        std::cout << "Node attributes: [" << element[0]["cell_set_designation"] << element[0]["taxonomy_id"] << "]" << std::endl;

        TaxonomyNode node("node_attributes");

        TaxonomyAttributes attributes;
        //attributes._map["members"]                 = element[0]["members"];
        //attributes._map["midpoint"]                = element[0]["midpoint"];
        //attributes._map["height"]                  = element[0]["height"];
        attributes._map["label"]                   = element[0]["label"];
        attributes._map["edgePar.col"]             = element[0]["edgePar.col"];
        //attributes._map["edgePar.lwd"]             = element[0]["edgePar.lwd"];
        //attributes._map["edgePar.conf"]            = element[0]["edgePar.conf"];
        attributes._map["cell_set_accession"]      = element[0]["cell_set_accession"];
        attributes._map["original_label"]          = element[0]["original_label"];
        attributes._map["cell_set_designation"]    = element[0]["cell_set_designation"];
        attributes._map["cell_set_alias"]          = element[0]["cell_set_alias"];
        attributes._map["cell_set_alt_alias"]      = element[0]["cell_set_alt_alias"];
        attributes._map["taxonomy_id"]             = element[0]["taxonomy_id"];
        attributes._map["comment"]                 = element[0]["comment"];
        //attributes._map["order"]                   = element[0]["order"];
        attributes._map["_row"]                    = element[0]["_row"];

        node.attributeIndex = taxonomy._nodeAttributes.size();
        taxonomy._nodeAttributes.push_back(attributes);

        return node;
    }

    TaxonomyNode parseLeafAttributes(Taxonomy& taxonomy, json element)
    {
        Q_ASSERT(element.size() == 1);

        std::cout << "Leaf attributes: [" << element[0]["cell_set_designation"] << element[0]["nodePar.col"] << "]" << std::endl;

        std::string name = element[0]["cell_set_designation"];
        TaxonomyNode node("leaf_attributes" + name);
        node.isLeaf = true;

        TaxonomyAttributes attributes;
        //attributes._map["members"]                 = element[0]["members"];
        //attributes._map["height"]                  = element[0]["height"];
        attributes._map["label"]                   = element[0]["label"];
        //attributes._map["leaf"]                    = element[0]["leaf"];
        attributes._map["edgePar.col"]             = element[0]["edgePar.col"];
        //attributes._map["edgePar.lwd"]             = element[0]["edgePar.lwd"];
        //attributes._map["edgePar.conf"]            = element[0]["edgePar.conf"];
        //attributes._map["nodePar.lab.cex"]         = element[0]["nodePar.lab.cex"];
        //attributes._map["nodePar.pch"]             = element[0]["nodePar.pch"];
        attributes._map["nodePar.lab.col"]         = element[0]["nodePar.lab.col"];
        //attributes._map["nodePar.cex"]             = element[0]["nodePar.cex"];
        attributes._map["nodePar.col"]             = element[0]["nodePar.col"];
        attributes._map["cell_set_accession"]      = element[0]["cell_set_accession"];
        attributes._map["original_label"]          = element[0]["original_label"];
        attributes._map["cell_set_designation"]    = element[0]["cell_set_designation"];
        attributes._map["cell_set_alias"]          = element[0]["cell_set_alias"];
        attributes._map["cell_set_alt_alias"]      = element[0]["cell_set_alt_alias"];
        attributes._map["taxonomy_id"]             = element[0]["taxonomy_id"];
        attributes._map["comment"]                 = element[0]["comment"];
        //attributes._map["order"]                   = element[0]["order"];
        attributes._map["_row"]                    = element[0]["_row"];

        node.attributeIndex = taxonomy._leafAttributes.size();
        taxonomy._leafAttributes.push_back(attributes);

        return node;
    }

    std::vector<TaxonomyNode> traverseJsonArray(Taxonomy& taxonomy, json element)
    {
        Q_ASSERT(element.is_array());

        qDebug() << "Size of array:" << element.size();

        std::vector<TaxonomyNode> arrayNodes;

        for (int i = 0; i < element.size(); i++)
        {
            qDebug() << "Children of element: " << i;
            TaxonomyNode aNode("array_node");
            std::vector<TaxonomyNode> nodes = traverseJsonObject(taxonomy, element[i]);

            for (TaxonomyNode& node : nodes)
                aNode._children.push_back(node);

            arrayNodes.push_back(aNode);
        }
        return arrayNodes;
    }

    std::vector<TaxonomyNode> traverseJsonObject(Taxonomy& taxonomy, json element)
    {
        Q_ASSERT(element.is_object());

        std::vector<TaxonomyNode> nodes;

        for (auto& el : element.items())
        {
            std::cout << "key: " << el.key() << std::endl;// << ", value:" << el.value() << '\n';

            TaxonomyNode parent(el.key());
            if (el.key() == "node_attributes")
            {
                nodes.push_back(parseNodeAttributes(taxonomy, el.value()));
                continue;
            }
            else if (el.key() == "leaf_attributes")
            {
                TaxonomyNode leaf = parseLeafAttributes(taxonomy, el.value());
                nodes.push_back(leaf);
                continue;
            }
            else if (el.value().is_array())
            {
                std::vector<TaxonomyNode> nodes = traverseJsonArray(taxonomy, el.value());

                for (TaxonomyNode& node: nodes)
                    parent._children.push_back(node);
            }
            else if (el.value().is_object())
            {
                std::vector<TaxonomyNode> nodes = traverseJsonObject(taxonomy, el.value());

                for (TaxonomyNode& node : nodes)
                    parent._children.push_back(node);
            }

            nodes.push_back(parent);
        }

        return nodes;
    }

    void printTreeNodeTraverse(TaxonomyNode& node, int level = 0)
    {
        for (int i = 0; i < level; i++)
            std::cout << "  ";
        std::cout << node._name << std::endl;
        for (TaxonomyNode& child : node._children)
        {
            printTreeNodeTraverse(child, level + 1);
        }
    }
}

Taxonomy Taxonomy::fromJsonFile()
{
    std::ifstream file("D:/Dropbox/Julian/Patchseq/taxonomy_2020.json");
    json data = json::parse(file);

    Taxonomy taxonomy;
    taxonomy._nodes = traverseJsonObject(taxonomy, data);

    taxonomy.printTree();

    return taxonomy;
}

void Taxonomy::printTree()
{
    for (TaxonomyNode& node: _nodes)
        printTreeNodeTraverse(node);
}

TaxonomyAttributes* Taxonomy::findLeafWithAttribute(std::string attributeName, std::string attributeValue)
{
    for (TaxonomyAttributes& attr : _leafAttributes)
    {
        if (attr._map.find(attributeName) != attr._map.end())
        {
            if (attr._map[attributeName] == attributeValue)
                return &attr;
        }
    }
    return nullptr;
}
