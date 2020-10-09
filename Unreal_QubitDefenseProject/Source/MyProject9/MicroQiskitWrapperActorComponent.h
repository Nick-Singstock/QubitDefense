// Fill out your copyright notice in the Description page of Project Settings.

#pragma once
#include <complex>
#include <cmath>
#include <vector>

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Math/IntPoint.h"
#include "MicroQiskitWrapperActorComponent.generated.h"

class UnrealQC {
public:
	int QubitDepth;
	int QubitCount;
	std::vector<std::complex<float>> QubitKetState;
	int mypow(int x, int p){
		if (p == 0) return 1;
		if (p == 1) return x;

		int tmp = mypow(x, p / 2);
		if ((p % 2) == 0) {
			return (tmp * tmp);
		}
		else {
			return (x * tmp * tmp);
		}
	}
	UnrealQC() {
		//QubitCount = qubits;
		//QubitDepth = std::complex<float>(0, 0);
		//
		//int LocalMaxIterate = mypow(2, QubitCount - 1);
		//for (int i = 0; i < LocalMaxIterate; i++) {
		//	QubitKetState[i] = std::complex<float>(0, 0);
		//}
		//QubitKetState[0] = std::complex<float>(1, 0);
	}
	int initialize(int qubits, int depth) {
		QubitCount = qubits;
		QubitDepth = 0;
		int LocalMaxIterate = mypow(2, QubitCount - 1);
		for (int i = 0; i < LocalMaxIterate; i++) {
			QubitKetState.push_back(std::complex<float>(0, 0));
		}
		QubitKetState[0] = std::complex<float>(1, 0);
		return 0;
	}

	int h(int m) {
		int LocalMaxIterate = mypow(2, QubitCount - 1);
		std::complex<float> s = std::complex<float>(0.70710678118, 0);
		for (int i = 0; i < QubitKetState.size(); i++) {
			int I = 2 * i - i % (mypow(2, m));
			int J = I + mypow(2, m);
			std::complex<float> a = s * QubitKetState[I] + s * QubitKetState[J];
			std::complex<float> b = s * QubitKetState[I] - s * QubitKetState[J];
			QubitKetState[I] = a;
			QubitKetState[J] = b;
		}
		return 0;
	}

	int x(int m) {
		int LocalMaxIterate = mypow(2, QubitCount - 1);
		for (int i = 0; i < QubitKetState.size(); i++) {
			int I = 2 * i - i % (mypow(2, m));
			int J = I + mypow(2, m);
			std::complex<float> a = QubitKetState[I];
			QubitKetState[I] = QubitKetState[J];
			QubitKetState[J] = a;
		}
		return 0;
	}


	int y(int m) {
		int LocalMaxIterate = mypow(2, QubitCount - 1);
		for (int i = 0; i < QubitKetState.size(); i++) {
			int I = 2 * i - i % (mypow(2, m));
			int J = I + mypow(2, m);
			std::complex<float> a = std::complex<float>(0, -1) * QubitKetState[I];
			QubitKetState[I] = std::complex<float>(0, 1) * QubitKetState[J];
			QubitKetState[J] = a;
		}
		return 0;
	}

	int z(int m) {
		int LocalMaxIterate = mypow(2, QubitCount - 1);
		for (int i = 0; i < QubitKetState.size(); i++) {
			int J = 2 * i - i % (mypow(2, m)) + mypow(2, m);
			QubitKetState[J] *= -1;
		}
		return 0;
	}
	int cnot(int c, int t) {

		int LocalMaxIterate = mypow(2, QubitCount - 2);
		for (int i = 0; i < LocalMaxIterate; i++) {
			int I = (mypow(2, c) + i % mypow(2, c) + ((i - i % mypow(2, c)) * 2) % mypow(2, t) + 2 * ((i - i % mypow(2, c)) * 2 - ((2 * (i - i % mypow(2, c))) % mypow(2, t))));
			int J = I + mypow(2, t);
			std::complex<float> Aux;
			Aux = QubitKetState[I];
			QubitKetState[I] = QubitKetState[J];
			QubitKetState[J] = Aux;
		}
		return 0;
	}

};

UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class MYPROJECT9_API UMicroQiskitWrapperActorComponent : public UActorComponent
{
	GENERATED_BODY()

public:	
	// Sets default values for this component's properties
	UMicroQiskitWrapperActorComponent();
protected:
	// Called when the game starts
	UnrealQC CoreUnrealQC;
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;
	//UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments")
	//	TMap<FIntPoint, FVector2D> GetStateKet(int InputOperationIndex);
	//UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<FVector2D> GetStateKetUnreal(/*TArray<FString> InputOperationIndex, TArray<int32> OperationIndex, TArray<int32> SecondOperationIndex*/);
	//UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<FVector2D> AddGate(FString GateType, int32 Input0, int32 input1);
	UPROPERTY(BlueprintReadOnly, VisibleAnywhere) TMap<FIntPoint, FVector2D> TMapRealImaginary;
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<FVector2D> CallQuantumHGate(int Input0);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<FVector2D> CallQuantumXGate(int Input0);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<FVector2D> CallQuantumYGate(int Input0);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<FVector2D> CallQuantumZGate(int Input0);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<FVector2D> CallQuantumCNOTGate(int Input0, int Input1);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") int CallQuantum(int Qubits);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<int> x_permut(int n, int m);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<int> cx_permut(int n, int c, int t);
};
