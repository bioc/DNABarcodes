// genericchromosome.h

#ifndef __GENERICCHROMOSOME_H_INCLUDED__ 
#define __GENERICCHROMOSOME_H_INCLUDED__

#include <ostream>

class GenericChromosome { 

  public:
    virtual               ~GenericChromosome();
    virtual unsigned int  fitness();
    virtual void          initialize();
    virtual unsigned int  getPrecalculatedFitness() const;
    virtual void          mutate();
    virtual bool          operator<(const GenericChromosome&) const;
    virtual bool          isFitnessKnown() const;
    virtual void          print(std::ostream& os) const;
    virtual GenericChromosome *clone() const;

    friend std::ostream& operator<<(std::ostream &os, GenericChromosome const& gc) {
      gc.print(os);
      return(os);
    }

  private:
    bool                  is_fitness_known_;
    unsigned int          fitness_;

};

#endif 
