#pragma once

class item{
public:
    uint32_t l;
    uint32_t m;
    uint32_t r;

    bool operator == (const item& element) const{
        return (this->l == element.l) && (this->m == element.m) && (this->r == element.r);
    }
};