// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Math/IntPoint.h"
#include "MicroQiskitWrapperActorComponent.generated.h"


UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class MYPROJECT9_API UMicroQiskitWrapperActorComponent : public UActorComponent
{
	GENERATED_BODY()

public:	
	// Sets default values for this component's properties
	UMicroQiskitWrapperActorComponent();

protected:
	// Called when the game starts
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;
	//UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments")
	//	TMap<FIntPoint, FVector2D> GetStateKet(int InputOperationIndex);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments")
		TArray<FVector2D> GetStateKet(int InputOperationIndex);

};
