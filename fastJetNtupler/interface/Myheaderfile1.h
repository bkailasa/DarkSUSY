#include "fastjet/PseudoJet.hh"



class Myheaderfile1 : public fastjet::PseudoJet::UserInfoBase{
 public:
   // default ctor
   //  - pdg_id        the PDG id of the particle
   //  - vertex_number theid of the vertex it originates from
   Myheaderfile1(const int & pdg_id_in, const int & vertex_number_in) : _pdg_id(pdg_id_in), _vertex_number(vertex_number_in){}
 
   /// access to the PDG id
   int pdg_id() const { return _pdg_id;}
   
   /// access to the vertex number
   int vertex_number() const { return _vertex_number;}
   
 protected:
 int _pdg_id;         // the associated pdg id
 int _vertex_number;  // the associated vertex number
 };
