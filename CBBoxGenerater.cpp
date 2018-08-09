#include "CBBoxGenerater.h"

const double CBBoxGenerater::TOL = 1.0e-10;//1,0*10^10

//�����U�s������
D3DXMATRIX CBBoxGenerater::GetCovarianceMatrix(
	D3DXVECTOR3* vertex,	// ���_�f�[�^�z�� 
	DWORD* index,		// �{�[���ɉe�����󂯂钸�_�C���f�b�N�X�̔z��
	int numVertices){	// �{�[���ɉe�����󂯂钸�_��

	float mx=0, my=0, mz=0;
	D3DXMATRIX retmat(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	//�{�[�����\�����钸�_�ɂ��ĕ��ϒl�����߂�
	for(int j=0; j<numVertices; j++){
		mx += vertex[index[j]].x;		// �w�l���v
		my += vertex[index[j]].y;		// �x�l���v
		mz += vertex[index[j]].z;		// �y�l���v
	}
	mx /= numVertices;		// �w�l����
	my /= numVertices;		// �x�l����
	mz /= numVertices;		// �y�l����
	//�����U�s��̊e���������߂�
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

//�ŗL�x�N�g�������߂�
bool CBBoxGenerater::GetEigenVector(D3DXMATRIX& a, D3DXMATRIX& x){
	int i,j,k,m,count;
	double amax, amax0, theta, co, si, co2, si2, cosi, pi = 4.0 * atan(1.0);
	double aii, aij, ajj, aik, ajk;
	bool result = false;
	// �P�ʍs��ɂ���
	D3DXMatrixIdentity(&x);
	count = 0;
	while(count <= MAX){
		//�Ίp�v�f�̍ő�l��T��
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
		//��������
		if(amax <= TOL){
			result = true;
			break;
		}else{
			aii = a.m[i][i]; 
			aij = a.m[i][j]; 
			ajj = a.m[j][j];
			//��]�p�x�v�Z
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

			//�����ϊ��s��
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
			//�ŗL�x�N�g��
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

//�ŗL�x�N�g������ɂ��Ăa�a�n�w���쐬
tag_BBOX CBBoxGenerater::CaclBBox(D3DXVECTOR3 *vertex, DWORD *index, int numVertices, D3DXMATRIX eigenmat){
	D3DXVECTOR3 rvec,svec,tvec,vertexvec,center;
	float rmin=0.0,rmax=0.0,smin=0.0,smax=0.0,tmin=0.0,tmax=0.0;
	float dot;
	tag_BBOX result;

	//�ŗL�x�N�g�����i�[���ꂽ�s����e�x�N�g���֑�������K��
	rvec.x = eigenmat._11; rvec.y = eigenmat._12; rvec.z = eigenmat._13;
	svec.x = eigenmat._21; svec.y = eigenmat._22; svec.z = eigenmat._23;
	tvec.x = eigenmat._31; tvec.y = eigenmat._32; tvec.z = eigenmat._33;
	//
	D3DXVec3Normalize(&rvec,&rvec);
	D3DXVec3Normalize(&svec,&svec);
	D3DXVec3Normalize(&tvec,&tvec);
	//�e�x�N�g���Ƃ��ׂĂ̒��_�̓��ς���ŏ��l�ƍő�l�����߂�
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
	//BBOX�\���̂ɒl���i�[
	result.vLocalPivot.x = rvec.x * center.x + svec.x * center.y + tvec.x * center.z;
	result.vLocalPivot.y = rvec.y * center.x + svec.y * center.y + tvec.y * center.z;
	result.vLocalPivot.z = rvec.z * center.x + svec.z * center.y + tvec.z * center.z;
	result.vLength.x = (rmax - rmin) / 2; result.vLength.y = (smax - smin) / 2; result.vLength.z = (tmax - tmin) / 2;
	result.vMax.x = rmax; result.vMax.y = smax; result.vMax.z = tmax;
	result.vMin.x = rmin; result.vMin.y = smin; result.vMin.z = tmin;
	result.vAxisX = rvec; result.vAxisY = svec; result.vAxisZ = tvec;
	return result;
}


