#include <vector>
#include <algorithm>

#include "AliPicoTrack.h"

bool sortPt(AliPicoTrack* t1, AliPicoTrack* t2) { return (t1->Pt() < t2->Pt()); }
bool sortEta(AliPicoTrack* t1, AliPicoTrack* t2) { return (t1->Eta() < t2->Eta()); }

void testPico()
{

  std::vector<AliPicoTrack*> vecT;
  AliPicoTrack* p1 = new AliPicoTrack(3,1.,1,2,0,0);
  AliPicoTrack* p2 = new AliPicoTrack(1.1,-2,0.0,0,0,0);
  AliPicoTrack* p3 = new AliPicoTrack(-2.1,1,0.0,0,0,0);
  AliPicoTrack* p4 = new AliPicoTrack(-1.1,2,0.0,0,0,0);

  vecT.push_back(p1);
  vecT.push_back(p2);
  vecT.push_back(p3);
  vecT.push_back(p4);

  printf("pre-sorting\n");
  for(int i = 0; i < vecT.size(); ++i) {
    printf("%g | %g\n",vecT[i]->Pt(),vecT[i]->Eta());
  }

  // // std::sort(vecT.begin(), vecT.end());
  std::sort(vecT.begin(), vecT.end(),sortEta);
  // std::sort(vecT.begin(), vecT.end(),sortEta);

  printf("post-sorting\n");
  for(int i = 0; i < vecT.size(); ++i) {
    printf("%g | %g\n",vecT[i]->Pt(),vecT[i]->Eta());
  }

};
