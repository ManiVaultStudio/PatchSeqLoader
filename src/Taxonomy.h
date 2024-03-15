#pragma once

#include <vector>
#include <unordered_map>
#include <string>

class TaxonomyAttributes
{
public:

public:
    std::unordered_map<std::string, std::string> _map;
};

class TaxonomyNode
{
public:
    TaxonomyNode() :
        _name("")
    {

    }

    TaxonomyNode(std::string name) :
        _name(name)
    {

    }

public:
    std::string _name;

    std::vector<TaxonomyNode> _children;

    int attributeIndex = -1;

    bool isLeaf = false;
};

class Taxonomy
{
public:
    static Taxonomy fromJsonFile();

    void printTree();

    TaxonomyAttributes* findLeafWithAttribute(std::string attributeName, std::string attributeValue);

    std::vector<TaxonomyNode> _nodes;

    std::vector<TaxonomyAttributes> _nodeAttributes;
    std::vector<TaxonomyAttributes> _leafAttributes;
};
