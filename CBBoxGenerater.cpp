#include "CBBoxGenerater.h"

const double CBBoxGenerater::TOL = 1.0e-10;//1,0*10^10

//共分散行列を作る
D3DXMATRIX CBBoxGenerater::GetCovarianceMatrix(
	D3DXVECTOR3* vertex,	// 頂点データ配列 
	DWORD* index,		// ボーンに影響を受ける頂点インデックスの配列
	int numVertices){	// ボーンに影響を受ける頂点数

	float mx=0, my=0, mz=0;
	D3DXMATRIX retmat(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	//ボーンを構成する頂点について平均値を求める
	for(int j=0; j<numVertices; j++){
		mx += vertex[index[j]].x;		// Ｘ値合計
		my += vertex[index[j]].y;		// Ｙ値合計
		mz += vertex[index[j]].z;		// Ｚ値合計
	}
	mx /= numVertices;		// Ｘ値平均
	my /= numVertices;		// Ｙ値平均
	mz /= numVertices;		// Ｚ値平均
	//共分散行列の各成分を求める
	for(int j=0; j<numVertices; j++){
		retmat._11 += (vertex[index[j]].x - mx) * (vertex[index[j]].x - mx);
		retmat._22 += (vertex[index[j]].y - my) * (vertex[index[j]].y - my);
		retmat._33 += (vertex[index[j]].z - mz) * (vertex[index[j]].z - mz);
		retmat._12 += (vertex[index[j]].x - mx) * (vertex[index[j]].y - my);
		retmat._13 += (vertex[index[j]].x - mx) * (vertex[index[j]].z - mz);
		retmat._23 += (vertex[index[j]].y - my) * (vertex[index[j]].z - mz);
	}
	retmat._11 /= numVertices;
	retmat._22 /= numVertices;
	retmat._33 /= numVertices;
	retmat._12 /= numVertices;
	retmat._13 /= numVertices;
	retmat._23 /= numVertices;
	retmat._21 = retmat._12;
	retmat._31 = retmat._13;
	retmat._32 = retmat._23;
	return retmat;
}

//固有ベクトルを求める
bool CBBoxGenerater::GetEigenVector(D3DXMATRIX& a, D3DXMATRIX& x){
	int i,j,k,m,count;
	double amax, amax0, theta, co, si, co2, si2, cosi, pi = 4.0 * atan(1.0);
	double aii, aij, ajj, aik, ajk;
	bool result = false;
	// 単位行列にする
	D3DXMatrixIdentity(&x);
	count = 0;
	while(count <= MAX){
		//対角要素の最大値を探索
		amax = 0.0;
		for(k=0; k<N-1; k++){
			for(m=k+1; m<N; m++){
				amax0 = fabs(*(a.m[N*k+m]));
				if(amax0 > amax){
					i = k;
					j = m;
					amax = amax0;
				}
			}
		}
		//収束判定
		if(amax <= TOL){
			result = true;
			break;
		}else{
			aii = a.m[i][i]; 
			aij = a.m[i][j]; 
			ajj = a.m[j][j];
			//回転角度計算
			if(fabs(aii-ajj) < TOL){
				theta  = 0.25 * pi * aij / fabs(aij);
			}
			else{
				theta = 0.5 * atan(2.0 * aij / (aii - ajj));
			}
			co = cos(theta); 
			si = sin(theta); 
			co2 = co * co; 
			si2 = si * si; 
			cosi = co * si;

			//相似変換行列
			a.m[i][i] = co2 * aii + 2.0 * cosi * aij + si2 * ajj;
			a.m[j][j] = si2 * aii - 2.0 * cosi * aij + co2 * ajj;
			a.m[i][j] = 0.0; 
			a.m[j][i] = 0.0;

			for(k=0; k<N; k++){
				if(k != i && k != j){
					aik = a.m[k][i];
					ajk = a.m[k][j];
					a.m[k][i] = co * aik + si * ajk;
					a.m[i][k] = a.m[k][i];
					a.m[k][j] = -si * aik + co * ajk;
					a.m[j][k] = a.m[k][j];
				}
			}
			//固有ベクトル
			for(k=0; k<N; k++){
				aik = x.m[k][i];
				ajk = x.m[k][j];
				x.m[k][i] = co * aik + si* ajk;
				x.m[k][j] = -si * aik + co * ajk;
			}
			count++;
		}
	}
	return result;
}

//固有ベクトルを基にしてＢＢＯＸを作成
tag_BBOX CBBoxGenerater::CaclBBox(D3DXVECTOR3 *vertex, DWORD *index, int numVertices, D3DXMATRIX eigenmat){
	D3DXVECTOR3 rvec,svec,tvec,vertexvec,center;
	float rmin=0.0,rmax=0.0,smin=0.0,smax=0.0,tmin=0.0,tmax=0.0;
	float dot;
	tag_BBOX result;

	//固有ベクトルが格納された行列を各ベクトルへ代入し正規化
	rvec.x = eigenmat._11; rvec.y = eigenmat._12; rvec.z = eigenmat._13;
	svec.x = eigenmat._21; svec.y = eigenmat._22; svec.z = eigenmat._23;
	tvec.x = eigenmat._31; tvec.y = eigenmat._32; tvec.z = eigenmat._33;
	//
	D3DXVec3Normalize(&rvec,&rvec);
	D3DXVec3Normalize(&svec,&svec);
	D3DXVec3Normalize(&tvec,&tvec);
	//各ベクトルとすべての頂点の内積から最小値と最大値を求める
	for(int i=0; i<numVertices; i++){
		vertexvec.x = vertex[index[i]].x; vertexvec.y = vertex[index[i]].y; vertexvec.z = vertex[index[i]].z;
		dot = D3DXVec3Dot(&rvec,&vertexvec);
		if(i==0 || dot < rmin) rmin = dot;
		if(i==0 || dot > rmax) rmax = dot;

		dot= D3DXVec3Dot(&svec,&vertexvec);
		if(i==0 || dot < smin) smin = dot;
		if(i==0 || dot > smax) smax = dot;

		dot= D3DXVec3Dot(&tvec,&vertexvec);
		if(i==0 || dot < tmin) tmin = dot;	
		if(i==0 || dot > tmax) tmax = dot;
	}
	center.x = (rmax + rmin) / 2;	center.y = (smax + smin) / 2;	center.z = (tmax + tmin) / 2;
	//BBOX構造体に値を格納
	result.vLocalPivot.x = rvec.x * center.x + svec.x * center.y + tvec.x * center.z;
	result.vLocalPivot.y = rvec.y * center.x + svec.y * center.y + tvec.y * center.z;
	result.vLocalPivot.z = rvec.z * center.x + svec.z * center.y + tvec.z * center.z;
	result.vLength.x = (rmax - rmin) / 2; result.vLength.y = (smax - smin) / 2; result.vLength.z = (tmax - tmin) / 2;
	result.vMax.x = rmax; result.vMax.y = smax; result.vMax.z = tmax;
	result.vMin.x = rmin; result.vMin.y = smin; result.vMin.z = tmin;
	result.vAxisX = rvec; result.vAxisY = svec; result.vAxisZ = tvec;
	return result;
}


