Int_t getTarget(Int_t fSystInd) {
  switch(fSystInd) {
    //Track cuts (dif. FB)
    case 1:
      return 0;
    case 9:
      return 1;
    case 10:
      return 1;
    case 11:
      return 1;
    case 12:
      return 1;
    case 13:
      return 0;//1;
    case 14:
      return 0;//1;
    case 15:
      return 0;//1;
    case 16:
      return 0;//1;
    //TPC nclusters
    case 6:
      return 1;//2;
    case 7:
      return 1;//2;
    case 8:
      return 1;//2;
    //DCAxy cuts -- disabled (due to eff)
    case 2:
      return 0;
    case 3:
      return 0;
    //DCAz cut:
    case 4:
      return 1;//2;//3;
    case 5:
      return 1;//2;//3;
    //vtx z position
    case 17:
      return 2;//3;//4;
    case 18:
      return 2;//3;//4;
    case 19:
      return 2;//3;//4;
    case 20:
      return 0;
    case 21:
      return 0;
    case 22:
      return 3;//4;//5;
    case 25:
      return 4;//5;//6;
    case 26:
      return 4;//5;//6;
    default:
      return 0;
  }
}
Int_t getCorrelationError(Int_t fSystInd) {
  switch(fSystInd) {
    //Track cuts (dif. FB)
    case 1:
      return 1;
    case 9:
      return 1;
    case 10:
      return 1;
    case 11:
      return 1;
    case 12:
      return 1;
    case 13:
      return 1;
    case 14:
      return 1;
    case 15:
      return 1;
    case 16:
      return 1;
    //DCAxy cuts -- disabled (due to eff)
    case 2:
      return 0;
    case 3:
      return 0;
    //DCAz cut:
    case 4:
      return -1;
    case 5:
      return -1;
    //TPC nclusters
    case 6:
      return -1;
    case 7:
      return -1;
    case 8:
      return -1;
    //vtx z position
    case 17:
      return -1;
    case 18:
      return -1;
    case 19:
      return -1;
    case 20:
      return -1;
    case 21:
      return -1;
    case 22:
      return -1;
    case 25:
      return -1;
    case 26:
      return -1;
    default:
      return 0;
  }
}
