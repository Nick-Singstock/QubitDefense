// Fill out your copyright notice in the Description page of Project Settings.


#include "MicroQiskitWrapperActorComponent.h"

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
TArray<FVector2D> UMicroQiskitWrapperActorComponent::GetStateKet(int InputOperationIndex)
{

	TArray<FVector2D> KetStateOutput;
	return KetStateOutput;
}