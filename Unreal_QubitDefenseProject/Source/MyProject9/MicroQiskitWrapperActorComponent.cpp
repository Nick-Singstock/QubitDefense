// Fill out your copyright notice in the Description page of Project Settings.


#include "MicroQiskitWrapperActorComponent.h"
//
//// Sets default values for this component's properties
//UMicroQiskitWrapperActorComponent::UMicroQiskitWrapperActorComponent()
//{
//	// Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
//	// off to improve performance if you don't need them.
//	PrimaryComponentTick.bCanEverTick = true;
//
//	// ...
//}
//
//
//// Called when the game starts
//void UMicroQiskitWrapperActorComponent::BeginPlay()
//{
//	Super::BeginPlay();
//
//	// ...
//	
//}
//
//
//// Called every frame
//void UMicroQiskitWrapperActorComponent::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
//{
//	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);
//
//	// ...
//}
////bool UMicroQiskitWrapperActorComponent::ReturnSuccess(TArray<FVector2D> InputSolution)
////{
////
////
////}
//
//TArray<int> UMicroQiskitWrapperActorComponent::x_permut(int n, int m) {
//    // performs permutation equivalent to X gate
//    // args:
//    //    n (int): number of qubits
//    //    m (int): qubit on which gate is applied
//    // rets:
//    //    perm (vector int): permutation of 2n elements
//
//    // initialize variables
//    TArray<int> perm;
//    for (int i = 0; i < (1 << n); i++) {
//        perm.Add(i);
//    }
//
//    // loop through entries
//    // avoid matrix multiplication, quite fancy footwork
//    for (int i = 0; i < (1 << (n - 1)); i++) {
//        int I = 2 * i - i % (1 << m);
//        int J = I + (1 << m);
//        int a = perm[I];
//        perm[I] = perm[J];
//        perm[J] = a;
//    }
//
//    return perm;
//}
//
//TArray<int> UMicroQiskitWrapperActorComponent::cx_permut(int n, int c, int t) {
//    // performs permutation equivalent to CX gate
//    // args:
//    //    n (int): number of qubits
//    //    c (int): control qubit
//    //    t (int): target qubit
//    // rets:
//    //    perm (vector int): permutation of 2n elements
//    TArray<int> perm;
//    int bitshiftN = (1 << n);
//    for (int i = 0; i < bitshiftN; i++) {
//        perm.Add(i);
//    }
//
//    // cf. x_permut
//    int bitshiftNM2 = (1 << (n - 2));
//    for (int i = 0; i < bitshiftNM2; i++) {
//        int I = (1 << c) + i % (1 << c) + ((i - i % (1 << c)) * 2) % (1 << t) + 2 * ((i - i % (1 << c)) * 2 - ((2 * (i - i % (1 << c)) % (1 << t))));
//        int J = I + (1 << t);
//        int a = perm[I];
//        perm[I] = perm[J];
//        perm[J] = a;
//    }
//
//    return perm;
//}

// Sets default values for this component's properties
UMicroQiskitWrapperActorComponent::UMicroQiskitWrapperActorComponent()
{
	// Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
	// off to improve performance if you don't need them.
	PrimaryComponentTick.bCanEverTick = true;

	// ...
}


// Called when the game starts
void UMicroQiskitWrapperActorComponent::BeginPlay()
{
	Super::BeginPlay();

	// ...

}


// Called every frame
void UMicroQiskitWrapperActorComponent::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	// ...
}
//bool UMicroQiskitWrapperActorComponent::ReturnSuccess(TArray<FVector2D> InputSolution)
//{
//
//
//}

TArray<int> UMicroQiskitWrapperActorComponent::x_permut(int n, int m) {
	// performs permutation equivalent to X gate
	// args:
	//    n (int): number of qubits
	//    m (int): qubit on which gate is applied
	// rets:
	//    perm (vector int): permutation of 2n elements

	// initialize variables
	TArray<int> perm;
	for (int i = 0; i < (1 << n); i++) {
		perm.Add(i);
	}

	// loop through entries
	// avoid matrix multiplication, quite fancy footwork
	for (int i = 0; i < (1 << (n - 1)); i++) {
		int I = 2 * i - i % (1 << m);
		int J = I + (1 << m);
		int a = perm[I];
		perm[I] = perm[J];
		perm[J] = a;
	}

	return perm;
}

TArray<int> UMicroQiskitWrapperActorComponent::cx_permut(int n, int c, int t) {
	// performs permutation equivalent to CX gate
	// args:
	//    n (int): number of qubits
	//    c (int): control qubit
	//    t (int): target qubit
	// rets:
	//    perm (vector int): permutation of 2n elements
	TArray<int> perm;
	int bitshiftN = (1 << n);
	for (int i = 0; i < bitshiftN; i++) {
		perm.Add(i);
	}

	// cf. x_permut
	int bitshiftNM2 = (1 << (n - 2));
	for (int i = 0; i < bitshiftNM2; i++) {
		int I = (1 << c) + i % (1 << c) + ((i - i % (1 << c)) * 2) % (1 << t) + 2 * ((i - i % (1 << c)) * 2 - ((2 * (i - i % (1 << c)) % (1 << t))));
		int J = I + (1 << t);
		int a = perm[I];
		perm[I] = perm[J];
		perm[J] = a;
	}

	return perm;
}



TArray<FVector2D> UMicroQiskitWrapperActorComponent::CallQuantumHGate(int Input0) {
	CoreUnrealQC.h(Input0);
	TArray<FVector2D> OutputFormat;
	for (int i = 0; i < CoreUnrealQC.QubitKetState.size(); i++) {
		FVector2D OutputValue = FVector2D(std::abs(CoreUnrealQC.QubitKetState[i]), std::arg(CoreUnrealQC.QubitKetState[i]));
		OutputFormat.Add(OutputValue);
	}
	return OutputFormat;
}
TArray<FVector2D> UMicroQiskitWrapperActorComponent::CallQuantumXGate(int Input0) {
	CoreUnrealQC.x(Input0);
	TArray<FVector2D> OutputFormat;
	for (int i = 0; i < CoreUnrealQC.QubitKetState.size(); i++) {
		FVector2D OutputValue = FVector2D(std::abs(CoreUnrealQC.QubitKetState[i]), std::arg(CoreUnrealQC.QubitKetState[i]));
		OutputFormat.Add(OutputValue);
}
return OutputFormat;
}
TArray<FVector2D> UMicroQiskitWrapperActorComponent::CallQuantumYGate(int Input0) {
	CoreUnrealQC.y(Input0);
	TArray<FVector2D> OutputFormat;
	for (int i = 0; i < CoreUnrealQC.QubitKetState.size(); i++) {
		FVector2D OutputValue = FVector2D(std::abs(CoreUnrealQC.QubitKetState[i]), std::arg(CoreUnrealQC.QubitKetState[i]));
		OutputFormat.Add(OutputValue);
	}
	return OutputFormat;
}

TArray<FVector2D> UMicroQiskitWrapperActorComponent::CallQuantumZGate(int Input0) {
	CoreUnrealQC.z(Input0);
	TArray<FVector2D> OutputFormat;
	for (int i = 0; i < CoreUnrealQC.QubitKetState.size(); i++) {
		FVector2D OutputValue = FVector2D(std::abs(CoreUnrealQC.QubitKetState[i]), std::arg(CoreUnrealQC.QubitKetState[i]));
		OutputFormat.Add(OutputValue);
	}
	return OutputFormat;
}

TArray<FVector2D> UMicroQiskitWrapperActorComponent::CallQuantumCNOTGate( int Input0, int Input1) {
	CoreUnrealQC.cnot(Input0, Input1);
	TArray<FVector2D> OutputFormat;
	for (int i = 0; i < CoreUnrealQC.QubitKetState.size(); i++) {
		FVector2D OutputValue = FVector2D(std::abs(CoreUnrealQC.QubitKetState[i]), std::arg(CoreUnrealQC.QubitKetState[i]));
		OutputFormat.Add(OutputValue);
	}
	return OutputFormat;
}

int UMicroQiskitWrapperActorComponent::CallQuantum(int QubitsInput) {
	CoreUnrealQC.initialize(QubitsInput, 0);
	return 0;
}