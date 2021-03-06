//----------------------------------------------------------------
// MargotDataRecordTree.cc
//
// - written by: Brian Roeder, LPC Caen, 15 Nov 08
// - email - broeder@comp.tamu.edu
//
// - Usage: A C++ class for GEANT4 for creating a ROOT tree for 
// -        the Margot neutron detector simulation event by event.
//
///////////////////////////////////////////////////////////////////
//
// 19 Jan 07 - For methods and functions associated with MargotDataRecordTree
//           - pointers, see MargotDataRecordTree.hh 
//
// 29 Jan 07 - Moved member methods and functions to this file "MargotDataRecordTree.cc"
//           - To keep better track of them
// 
//
#include "Randomize.hh" 
#include "MargotDataRecordTree.hh"

// Access to Analysis pointer! (see MargotSD.cc EndOfEvent() for Example)
MargotDataRecordTree* MargotDataRecordTree::MargotPointer;

MargotDataRecordTree::MargotDataRecordTree() : 
  // Initialized Values
  eng_int(0), theta_int(0), phi_int(0), 
  eng_MicroBulk_Tot(0), X_MicroBulk(0), Y_MicroBulk(0), Z_MicroBulk(0), PosMag_MicroBulk(0),
  eng_BulkFront(0), X_BulkFront(0), Y_BulkFront(0), Z_BulkFront(0), PosMag_BulkFront(0),
  eng_BulkBack(0), X_BulkBack(0), Y_BulkBack(0), Z_BulkBack(0), PosMag_BulkBack(0),
  engDepEvent(0), PosMagEvent(0),
  event_counter(0), number_detected(0), number_MicroBulk(0),
  number_BulkFront(0), number_BulkBack(0),number_PunchThrough(0)
{ /* Constructor */
 // When Pointer is constructed, assigns address of this class to it.
  MargotPointer = this;

  for(G4int i=0;i<5;i++)
    {eng_MicroBulk[i] = 0.;}

  // Create pointers to ROOT Analysis Tree and Branches

  cout <<"\n Starting Data Tree Constructor" << endl;
  
  const Char_t* evt_file = "IonChamberDataTree.root";

  DataFile = new TFile(evt_file,"RECREATE");
  
  IonChDistTree = new TTree("IonChDistTree","Ion Chamber Data");
  IonChDistTree->Branch("eng_MicroBulk1",&eng_MicroBulk[0],"eng_MicroBulk/D");
  IonChDistTree->Branch("eng_MicroBulk2",&eng_MicroBulk[1],"eng_MicroBulk/D");
  IonChDistTree->Branch("eng_MicroBulk3",&eng_MicroBulk[2],"eng_MicroBulk/D");
  IonChDistTree->Branch("eng_MicroBulk4",&eng_MicroBulk[3],"eng_MicroBulk/D");
  IonChDistTree->Branch("eng_MicroBulk5",&eng_MicroBulk[4],"eng_MicroBulk/D");
  IonChDistTree->Branch("eng_MicroBulk_Tot",&eng_MicroBulk_Tot,"eng_MicroBulk_Tot/D");
  IonChDistTree->Branch("X_MicroBulk",&X_MicroBulk,"X_MicroBulk/D");
  IonChDistTree->Branch("Y_MicroBulk",&Y_MicroBulk,"Y_MicroBulk/D");
  IonChDistTree->Branch("Z_MicroBulk",&Z_MicroBulk,"Z_MicroBulk/D");
  IonChDistTree->Branch("PosMag_MicroBulk",&PosMag_MicroBulk,"PosMag_MicroBulk/D");

  IonChDistTree->Branch("eng_BulkFront",&eng_BulkFront,"eng_BulkFront/D");
  IonChDistTree->Branch("X_BulkFront",&X_BulkFront,"X_BulkFront/D");
  IonChDistTree->Branch("Y_BulkFront",&Y_BulkFront,"Y_BulkFront/D");
  IonChDistTree->Branch("Z_BulkFront",&Z_BulkFront,"Z_BulkFront/D");
  IonChDistTree->Branch("PosMag_BulkFront",&PosMag_BulkFront,"PosMag_BulkFront/D");

  IonChDistTree->Branch("eng_BulkBack",&eng_BulkBack,"eng_BulkBack/D");
  IonChDistTree->Branch("X_BulkBack",&X_BulkBack,"X_BulkBack/D");
  IonChDistTree->Branch("Y_BulkBack",&Y_BulkBack,"Y_BulkBack/D");
  IonChDistTree->Branch("Z_BulkBack",&Z_BulkBack,"Z_BulkBack/D");
  IonChDistTree->Branch("PosMag_BulkBack",&PosMag_BulkBack,"PosMag_BulkBack/D");

  h1_ZStopPos = new TH1F("ZStopPos","Z_Stop Position in MicroBulk (mm)",241,0.,240.);
  h1_ZStopPos->GetXaxis()->SetTitle("mm, New Center=120.0mm");
  h1_ZStopPos->GetYaxis()->SetTitle("Counts/mm"); 

  cout << MargotPointer << endl;
}

MargotDataRecordTree::~MargotDataRecordTree()
{/* Destructor, Close root file */
 
  DataFile->Write(); 
  DataFile->Close();
  cout << "############################################################" << endl;
  cout << "Created Root Tree File = \"IonChamberDataTree.root\"" << endl;
  cout << "Got the ROOT Tree File from Data Record Routine!" << endl; 
  delete DataFile;
}

// Member functions of MargotDataRecordTree class

void MargotDataRecordTree::senddataPG(double value1=0.,G4double theta=0., G4double phi=0.)
{
    eng_int = value1;
    event_counter++;
    theta_int = theta/deg;
    phi_int = phi/deg;
}
 
void MargotDataRecordTree::ShowDataFromEvent()
{

      G4cout << "====================================================" << G4endl;
      G4cout << "The Energy of the Initial Particle was :    " << G4BestUnit(eng_int,"Energy") << G4endl;
      G4cout << "The Z - Position of the Hit was        :    " << PosMagEvent << " mm" << G4endl;
      G4cout << "Measured Eng Dep. in Detectors         :    " << G4BestUnit(engDepEvent, "Energy") << G4endl;
      G4cout << "Total Detected So Far                  :    " << number_detected << G4endl;
      G4cout << "Number Seen in BulkFront So Far        :    " << number_BulkFront << G4endl;
      G4cout << "Number Seen in MicroBulk So Far        :    " << number_MicroBulk << G4endl;
      G4cout << "Number Seen in BulkBack So Far         :    " << number_BulkBack << G4endl;
      G4cout << "Total Punch Through Detector So Far    :    " << number_PunchThrough << G4endl;
      G4cout << endl;     
}

void MargotDataRecordTree::GetParticleTotals()
{
  G4int number_stopped = 0;
  
  if(number_PunchThrough > 0)
    {number_stopped = number_detected - number_PunchThrough;}
  else if (number_PunchThrough <= 0)
    {number_stopped = number_detected;}


      G4cout << "The Initial Number of Beam Particles was           : " << event_counter << G4endl;
      G4cout << "Number of Beam Particles Detected                  : " << number_detected << G4endl;
      G4cout << "Number of Beam Particles Stopped in Detector Region: " << number_stopped << G4endl;
      G4cout << "Number of Beam Particles Seen in BulkFront         : " << number_BulkFront << G4endl;
      G4cout << "Number of Beam Particles Seen in MicroBulk         : " << number_MicroBulk << G4endl;
      G4cout << "Number of Beam Particles Seen in BulkBack          : " << number_BulkBack << G4endl;
      G4cout << "Number of Beam Particles Punch Through Det.        : " << number_PunchThrough << G4endl;
      G4cout << G4endl;
}
 
void MargotDataRecordTree::sendDeltaEData(G4double MicroBulk_DE[], G4ThreeVector Pos1,G4double MicroBulk_DE_Tot, G4double DE2, G4ThreeVector Pos2, G4double DE3, G4ThreeVector Pos3)
{
  // In IonChamberSD.cc
  // "MicroBulk_DE[]" is the energy loss in each MicroBulk pad
  // "MicroBulk_DE_Tot is the sum of all MicroPads (added separately).
  // "DE2" is the Bulk Front Det.
  // "DE3" is the Bulk Back Det. 

  if(Pos3(2) >= 170.0*mm)
    {number_PunchThrough++;}

  if(MicroBulk_DE_Tot > 100.*keV)
    {  
      number_detected++;
      // MicroBulk
      if(MicroBulk_DE_Tot > 100.*keV)
	{number_MicroBulk++;}
      eng_MicroBulk[0] = MicroBulk_DE[0];
      eng_MicroBulk[1] = MicroBulk_DE[1];
      eng_MicroBulk[2] = MicroBulk_DE[2];
      eng_MicroBulk[3] = MicroBulk_DE[3];
      eng_MicroBulk[4] = MicroBulk_DE[4];
      eng_MicroBulk_Tot = MicroBulk_DE_Tot;
      X_MicroBulk = Pos1(0);
      Y_MicroBulk = Pos1(1);
      Z_MicroBulk = Pos1(2);
  
      // Binning fixer
      G4double FixBinRand = (G4UniformRand()-0.5)*mm;

      if(Z_MicroBulk >= 70*mm && Z_MicroBulk < 170.0*mm)
	{h1_ZStopPos->Fill(Z_MicroBulk+FixBinRand);}
      else if(Pos2(2) < 70*mm)
	{h1_ZStopPos->Fill(Pos2(2)+FixBinRand);}
      else if(Pos3(2) >= 170.0*mm)
	{h1_ZStopPos->Fill(Pos3(2)+FixBinRand);}

      G4double x1_sq = pow(X_MicroBulk,2);
      G4double y1_sq = pow(Y_MicroBulk,2);
      G4double z1_sq = pow(Z_MicroBulk,2);

      PosMag_MicroBulk = sqrt(x1_sq+y1_sq+z1_sq);

      // Bulk Front
      if(DE2 > 100.*keV)
	{number_BulkFront++;}
      eng_BulkFront = DE2;
      X_BulkFront = Pos2(0);
      Y_BulkFront = Pos2(1);
      Z_BulkFront = Pos2(2);
  
      G4double x2_sq = pow(X_BulkFront,2);
      G4double y2_sq = pow(Y_BulkFront,2);
      G4double z2_sq = pow(Z_BulkFront,2);

      PosMag_BulkFront = sqrt(x2_sq+y2_sq+z2_sq);

      // Bulk Back
      if(DE3 > 100.*keV)
	{number_BulkBack++;}
      eng_BulkBack = DE3;
      X_BulkBack = Pos3(0);
      Y_BulkBack = Pos3(1);
      Z_BulkBack = Pos3(2);
  
      G4double x3_sq = pow(X_BulkBack,2);
      G4double y3_sq = pow(Y_BulkBack,2);
      G4double z3_sq = pow(Z_BulkBack,2);

      PosMag_BulkBack = sqrt(x3_sq+y3_sq+z3_sq);

      IonChDistTree->Fill();
    }

  // Event Data Holders (for display only)
  engDepEvent = MicroBulk_DE_Tot+DE2+DE3;
  if(DE2 <= 100.*keV)    
    {PosMagEvent = 0.;}
  if(DE2 > 100.*keV)          // DE 2 is Bulk Front
    {PosMagEvent = Pos2(2);}
  if(eng_MicroBulk_Tot > 100.*keV) // MicroBulk
    {PosMagEvent = Pos1(2);}
  if(DE3 > 100.*keV)          // DE3 is Bulk Back
    {PosMagEvent = Pos3(2);}
 
  // Zero the Data Holders for the Next Event

  for(G4int i=0;i<5;i++)
    {eng_MicroBulk[i] = 0.;}

  eng_MicroBulk_Tot = 0.;
  X_MicroBulk = 0.;
  Y_MicroBulk = 0.;
  Z_MicroBulk = 0.;

  PosMag_MicroBulk = 0.;

  eng_BulkFront = 0.;
  X_BulkFront = 0.;
  Y_BulkFront = 0.;
  Z_BulkFront = 0.;

  PosMag_BulkFront = 0.;

  eng_BulkBack = 0.;
  X_BulkBack = 0.;
  Y_BulkBack = 0.;
  Z_BulkBack = 0.;

  PosMag_BulkBack = 0.;

}

// End of DataRecordTree Functions!
